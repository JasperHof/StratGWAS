// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppParallel)]]
#include <RcppEigen.h>
#include <Rcpp.h>
#include <RcppParallel.h>
#include <Eigen/Eigenvalues>
#include "readBedBlock.h"
#include "geno_utils.h"
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <cmath>
#include <utility>
#include <cstdint>
#include <map>

using namespace Rcpp;

// ===========================================================================
// StratGWAS -- combined-panel chunk HE + saddlepoint, with KVIK LOCO offsets
// and a large-window heritability calibration.
//
// ONE genotype file holds every variant (WES/WGS rare AND array/imputed
// common). An annotation matrix assigns each SNP to exactly one CATEGORY
// (MAF bin, functional class, "common", ...). Every category is a variance
// component; the moment system is solved jointly across the categories present
// in each unit of analysis.
//
// TWO PASSES over the same file, same kernels, same solver:
//
//   PASS 1 -- ASSOCIATION.  Unit = a chunk of `chunk_size` consecutive SNPs.
//     No flanks (as SAIGE-GENE: a gene is tested marginally). Phenotype for
//     chromosome c is  y - LOCO_c , the LDAK-KVIK step-1 leave-one-chromosome-
//     out PRS, so the polygenic background of every OTHER chromosome is removed
//     while the tested chromosome's own signal is untouched. Output: vg, se_vg,
//     h2 and an SPA p-value per chunk x category x trait.
//
//   PASS 2 -- HERITABILITY.  Unit = a window of ~window_bp built from N whole
//     chunks, with the neighbouring windows as flanks. Phenotype = y - FULL PRS
//     (the full PRS is recovered exactly from the LOCO file, see prs_load).
//     Per category, the PRS actually subtracted is chosen by `prs_mask`, so a
//     rare-variant category can be estimated CONDITIONAL on the common-SNP PRS
//     while a common category is estimated on the raw phenotype. This is the
//     configuration that gave unbiased estimates in simulation; the chunk
//     estimates of pass 1 do not (LD leakage at 256-SNP resolution).
//
//   ADJUSTMENT.  Within each window, (window vg - sum of its chunks' vg) is
//   distributed back over the chunks in proportion to each chunk's expected
//   variance share under the alpha model, sum_j [2 f_j (1-f_j)]^(1+alpha).
//   p-values are NOT touched; vg_adj / h2_adj are reported next to vg / h2.
//
// Everything numerical below the context/driver layer (bed reader, cumulant
// algebra, quad_spa_solve, bern_spa_solve, the co-heritability machinery and
// test_chunk_annot's saddlepoint) is carried over verbatim from the validated
// single-file version. Removed: the separate common fileset, off_diag,
// project_common, per-SNP weights, flank_categories, the GRM background
// component, explicit `pairs`, binary_gap, and the no-annotation code path
// (an absent annotation now means one category, "ALL").
// ===========================================================================

namespace {


typedef Eigen::MatrixXf GenoMat;

struct BimInfo {
    std::vector<std::string> chr;
    std::vector<long>        bp;
    int n_snps;
};

static BimInfo read_bim_positions(const std::string& bim_path) {
    BimInfo out;
    std::ifstream in(bim_path.c_str());
    if (!in.is_open()) stop("Could not open .bim file: " + bim_path);
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::istringstream ss(line);
        std::string chr, rsid, cm, bp;
        ss >> chr >> rsid >> cm >> bp;
        out.chr.push_back(chr);
        out.bp.push_back(std::stol(bp));
    }
    out.n_snps = (int) out.bp.size();
    return out;
}

static int lower_index(const std::vector<long>& bp, int lo, int hi, long value) {
    std::vector<long>::const_iterator b = bp.begin();
    std::vector<long>::const_iterator it = std::lower_bound(b + lo, b + hi + 1, value);
    return (int)(it - b);
}

// ===========================================================================
// THREAD-SAFE .bed READER
//
// The previous design read genotypes through readBedBlock(), which returns an
// Rcpp IntegerMatrix. That has two costs, and the second one is what forced the
// whole assembly stage to be single-threaded:
//   * an n_total x m block of int is allocated and filled (467 MB per chunk at
//     n = 1e5, m = 1168) and then gathered into float through a strided index,
//     so every genotype is written twice and read three times; and
//   * it touches R's heap, which is NOT thread-safe, so make_chunk() could only
//     ever run on the main thread while every worker sat idle.
//
// This reader decodes the 2-bit PLINK representation straight into the
// destination float matrix, for the kept individuals only, using a private
// ifstream per call. No R API, no intermediate, no gather -> assembly can move
// into the parallel workers.
//
// DECODE TABLE. Rather than hard-code a genotype convention (which differs
// between readers, and a silent mismatch would flip alleles), the 4-entry table
// is CALIBRATED at start-up against readBedBlock() itself and then VERIFIED on
// random blocks elsewhere in the file. If they disagree the run stops. So this
// reader reproduces the old one bit-for-bit by construction.
// ===========================================================================
struct BedReader {
    std::string path;
    int  n_total;
    int  n_snps;
    long bytes_per_snp;
    float code2val[4];
    bool ready;
    BedReader() : n_total(0), n_snps(0), bytes_per_snp(0), ready(false) {
        // PLINK-1 defaults; overwritten by calibrate() before any real use.
        code2val[0] = 0.f; code2val[1] = -1.f; code2val[2] = 1.f; code2val[3] = 2.f;
    }
    void init(const std::string& prefix, int ntot, int nsnp) {
        path = prefix + ".bed"; n_total = ntot; n_snps = nsnp;
        bytes_per_snp = ((long) ntot + 3) / 4; ready = true;
    }
    // Raw (undecoded-to-float) read of SNPs [i0,i1) for individuals `keep`.
    // Returns false on any I/O problem -- callers must NOT stop() from a worker
    // thread, so failure is propagated as a skipped chunk instead.
    bool read_raw(const std::vector<int>& keep, int i0, int i1, GenoMat& X) const {
        const int m = i1 - i0, nk = (int) keep.size();
        X.resize(nk, m > 0 ? m : 0);
        if (m <= 0) return true;
        if (!ready || i0 < 0 || i1 > n_snps) return false;
        std::ifstream in(path.c_str(), std::ios::binary);
        if (!in.is_open()) return false;
        std::vector<unsigned char> buf((size_t) m * (size_t) bytes_per_snp);
        in.seekg((std::streamoff) 3 + (std::streamoff) i0 * bytes_per_snp, std::ios::beg);
        in.read((char*) buf.data(), (std::streamsize) buf.size());
        if (in.gcount() != (std::streamsize) buf.size()) return false;
        for (int j = 0; j < m; ++j) {
            const unsigned char* row = buf.data() + (size_t) j * (size_t) bytes_per_snp;
            float* col = X.data() + (size_t) j * (size_t) nk;   // Eigen is column-major
            for (int i = 0; i < nk; ++i) {
                const int gi = keep[i];
                col[i] = code2val[(row[gi >> 2] >> ((gi & 3) << 1)) & 3];
            }
        }
        return true;
    }
    // Decode an ARBITRARY (sorted) list of SNP indices. Used when a bp window has
    // to be sub-sampled: the caller picks the columns first, so the decoded float
    // matrix is bounded by the number of columns actually wanted rather than by
    // however many SNPs happen to fall in the window. One seek+read per SNP, each
    // bytes_per_snp long, so the byte buffer is bounded too.
    bool read_raw_sel(const std::vector<int>& keep, const std::vector<int>& snps,
                      GenoMat& X) const {
        const int m = (int) snps.size(), nk = (int) keep.size();
        X.resize(nk, m > 0 ? m : 0);
        if (m <= 0) return true;
        if (!ready) return false;
        std::ifstream in(path.c_str(), std::ios::binary);
        if (!in.is_open()) return false;
        std::vector<unsigned char> row((size_t) bytes_per_snp);
        for (int j = 0; j < m; ++j) {
            const int s = snps[j];
            if (s < 0 || s >= n_snps) return false;
            in.seekg((std::streamoff) 3 + (std::streamoff) s * bytes_per_snp, std::ios::beg);
            in.read((char*) row.data(), (std::streamsize) bytes_per_snp);
            if (in.gcount() != (std::streamsize) bytes_per_snp) return false;
            float* col = X.data() + (size_t) j * (size_t) nk;
            for (int i = 0; i < nk; ++i) {
                const int gi = keep[i];
                col[i] = code2val[(row[gi >> 2] >> ((gi & 3) << 1)) & 3];
            }
        }
        return true;
    }
};

// Learn the 4-entry decode table from readBedBlock(), then verify it on blocks
// scattered through the file. Main thread only (calls the R API).
static void calibrate_bed_reader(BedReader& br, const std::string& prefix,
                                 int n_total, int n_snps, const char* label) {
    br.init(prefix, n_total, n_snps);
    if (n_snps <= 0 || n_total <= 0) return;
    const int probe = std::min(n_snps, 64);
    std::vector<int> all(n_total);
    for (int i = 0; i < n_total; ++i) all[i] = i;

    // ---- calibrate on the first `probe` SNPs -------------------------------
    IntegerMatrix ref = readBedBlock(prefix + ".bed", n_total, n_snps, 0, n_total - 1, 0, probe - 1);
    std::ifstream in((prefix + ".bed").c_str(), std::ios::binary);
    if (!in.is_open()) stop("Could not open .bed file: " + prefix + ".bed");
    std::vector<unsigned char> buf((size_t) probe * (size_t) br.bytes_per_snp);
    in.seekg(3, std::ios::beg);
    in.read((char*) buf.data(), (std::streamsize) buf.size());
    if (in.gcount() != (std::streamsize) buf.size())
        stop("Short read while calibrating the .bed reader for " + prefix);

    bool seen[4] = {false,false,false,false};
    float val[4]  = {0.f,0.f,0.f,0.f};
    for (int j = 0; j < probe; ++j) {
        const unsigned char* row = buf.data() + (size_t) j * (size_t) br.bytes_per_snp;
        for (int i = 0; i < n_total; ++i) {
            int code = (row[i >> 2] >> ((i & 3) << 1)) & 3;
            float v  = (float) ref(i, j);
            if (!seen[code]) { seen[code] = true; val[code] = v; }
            else if (val[code] != v)
                stop("The .bed 2-bit code " + std::to_string(code) +
                     " maps to two different values in " + prefix +
                     " -- readBedBlock() is not a pure per-genotype decode, so the fast "
                     "reader cannot replace it. Please share readBedBlock.h.");
        }
    }
    for (int c = 0; c < 4; ++c) if (seen[c]) br.code2val[c] = val[c];

    // ---- verify on blocks elsewhere in the file ----------------------------
    int checks = 0;
    for (int rep = 1; rep <= 3 && n_snps > probe; ++rep) {
        long s = (long) n_snps * rep / 4;
        int i0 = (int) std::min<long>(s, (long) n_snps - probe);
        if (i0 <= 0) continue;
        IntegerMatrix r2 = readBedBlock(prefix + ".bed", n_total, n_snps, 0, n_total - 1, i0, i0 + probe - 1);
        GenoMat X;
        if (!br.read_raw(all, i0, i0 + probe, X))
            stop("Fast .bed reader failed a verification read on " + prefix);
        for (int j = 0; j < probe; ++j)
            for (int i = 0; i < n_total; ++i)
                if (X(i, j) != (float) r2(i, j))
                    stop("Fast .bed reader disagrees with readBedBlock() on " + prefix +
                         " -- refusing to run. Please share readBedBlock.h.");
        ++checks;
    }
    Rcout << "Fast .bed reader calibrated for " << label
          << " (codes 0/1/2/3 -> " << br.code2val[0] << "/" << br.code2val[1]
          << "/" << br.code2val[2] << "/" << br.code2val[3]
          << "; verified on " << checks << " extra blocks)\n";
}

// A read + standardized genotype cell: the standardized columns PLUS the per-
// column MAF (needed for alpha weighting) and the .bim index of its first SNP
// (needed to look up categories for the middle window).
struct Cell {
    GenoMat X;                     // n_keep x m, standardized (unweighted)
    std::vector<float> maf;        // per-column folded MAF (0 for collapsed pseudo-markers)
    int i0;                        // .bim index of column 0 (-1 if empty)
    std::vector<int> cat;          // per-column category id (only set when collapsing;
                                   // empty -> columns map 1:1 to .bim via i0)
};

// Standardize columns to unit variance (mean-impute missing == -1), capturing
// MAF. Reductions in double.
static void standardize_capture_maf(GenoMat& X, std::vector<float>& maf) {
    const int n = (int) X.rows();
    maf.assign(X.cols(), 0.f);
    for (int col = 0; col < X.cols(); ++col) {
        double sum = 0.0; int n_valid = 0;
        for (int i = 0; i < n; ++i) { float g = X(i, col); if (g != -1.f) { sum += g; ++n_valid; } }
        if (n_valid == 0) { X.col(col).setZero(); maf[col] = 0.f; continue; }
        double mean = sum / n_valid;
        double f = mean / 2.0; if (f > 0.5) f = 1.0 - f;   // fold to minor allele
        maf[col] = (float) f;
        double sq = 0.0;
        for (int i = 0; i < n; ++i) { double g = (X(i, col) == -1.f) ? mean : (double) X(i, col); sq += (g - mean) * (g - mean); }
        double sd = std::sqrt(sq / (n - 1));
        if (sd > 1e-10) { for (int i = 0; i < n; ++i) { double g = (X(i, col) == -1.f) ? mean : (double) X(i, col); X(i, col) = (float)((g - mean) / sd); } }
        else X.col(col).setZero();
    }
}

// Read a WES cell by SNP-INDEX range [i0, i1) (a cell is defined by SNP indices,
// not a bp grid, now that cells have variable width -- see build_cells).
// Thread-safe: returns false on I/O failure rather than calling stop().
static bool read_cell_idx(const BedReader& br, const std::vector<int>& keep,
                          int i0, int i1, Cell& cell,
                          const std::vector<int>* cat_id = 0,
                          int collapse_mac = 0, int collapse_n = 5) {
    cell = Cell(); cell.i0 = -1;
    int m = i1 - i0;
    int n_keep = (int) keep.size();
    if (m <= 0) { cell.X = GenoMat(n_keep, 0); return true; }
    GenoMat block;                                  // n_keep x m, RAW codes
    if (!br.read_raw(keep, i0, i1, block)) return false;

    // -------- no collapsing: original path (columns map 1:1 to the .bim) --------
    if (!(collapse_mac > 0 && cat_id != 0 && collapse_n >= 1)) {
        standardize_capture_maf(block, cell.maf);   // standardize in place
        cell.X = std::move(block);
        cell.i0 = i0;
        return true;
    }

    // -------- collapsing: variants with MAC < collapse_mac are grouped in runs --
    // of up to collapse_n CONSECUTIVE ultra-rare variants sharing the same category,
    // and additively collapsed (sum of MINOR-allele dosages) into one pseudo-marker.
    // Variants with MAC >= collapse_mac stay individual. Pseudo-markers carry MAF 0
    // so they are exempt from alpha weighting (their frequency is not well defined).
    int ncat_all = (int) cat_id->size();
    std::vector< std::vector<float> > out_cols;   // each length n_keep (raw scale)
    std::vector<int>  out_cat;                     // per output column: category id
    std::vector<char> out_pseudo;                  // per output column: 1 = collapsed

    std::vector<float> grp; int grp_cat = -1, grp_cnt = 0; bool grp_open = false;
    // flush the pending group into one pseudo-marker column
    // (declared as a lambda capturing the accumulators)
    // (kept inline below to avoid std::function overhead)

    for (int j = 0; j < m; ++j) {
        // fold to minor-allele dosage and compute MAC (missing excluded from MAC,
        // treated as 0 dosage for the burden)
        double s = 0.0; int nv = 0;
        for (int i = 0; i < n_keep; ++i) { float g = block(i, j); if (g >= 0.f) { s += g; ++nv; } }
        int si = (int) (s + 0.5);
        int mac = (nv > 0) ? std::min(si, 2 * nv - si) : 0;
        bool flip = (nv > 0 && s > (double) nv);          // coded allele is major
        int cid = (i0 + j >= 0 && i0 + j < ncat_all) ? (*cat_id)[i0 + j] : -1;

        if (mac >= collapse_mac) {
            if (grp_open) { out_cols.push_back(grp); out_cat.push_back(grp_cat); out_pseudo.push_back(1);
                            grp_open = false; grp_cnt = 0; grp_cat = -1; }
            std::vector<float> col(n_keep);
            for (int i = 0; i < n_keep; ++i) { float g = block(i, j); col[i] = (g < 0.f) ? -1.f : g; }
            out_cols.push_back(std::move(col)); out_cat.push_back(cid); out_pseudo.push_back(0);
        } else {
            if (grp_open && (grp_cat != cid || grp_cnt >= collapse_n)) {
                out_cols.push_back(grp); out_cat.push_back(grp_cat); out_pseudo.push_back(1);
                grp_open = false; grp_cnt = 0; grp_cat = -1;
            }
            if (!grp_open) { grp.assign(n_keep, 0.f); grp_open = true; grp_cat = cid; grp_cnt = 0; }
            for (int i = 0; i < n_keep; ++i) {
                float g = block(i, j);
                if (g >= 0.f) grp[i] += (flip ? (2.f - g) : g);   // additive minor dosage
            }
            ++grp_cnt;
        }
    }
    if (grp_open) { out_cols.push_back(grp); out_cat.push_back(grp_cat); out_pseudo.push_back(1); }

    int nout = (int) out_cols.size();
    GenoMat X(n_keep, nout);
    for (int k = 0; k < nout; ++k)
        for (int i = 0; i < n_keep; ++i) X(i, k) = out_cols[k][i];
    standardize_capture_maf(X, cell.maf);
    for (int k = 0; k < nout; ++k) if (out_pseudo[k]) cell.maf[k] = 0.f;   // exempt from alpha
    cell.X = std::move(X);
    cell.cat = out_cat;
    cell.i0 = i0;
    return true;
}

static void project_covariates(GenoMat& X, const GenoMat& Z, const GenoMat& M) {
    if (X.cols() == 0 || Z.cols() == 0) return;
    X.noalias() -= Z * (M * X);
    const int n = (int) X.rows();
    for (int j = 0; j < X.cols(); ++j) {
        double ss = X.col(j).cast<double>().squaredNorm();
        double sd = std::sqrt(ss / (n - 1));
        if (sd > 1e-10) X.col(j) /= (float) sd; else X.col(j).setZero();
    }
}

// Project a chunk's COMMON block out of a WES block:  X <- X - C (C'C)^+ C'X,
// then restore unit column variance (same convention as project_covariates).
// The WES kernels then cannot see any signal the local common SNPs can explain,
// so the estimand becomes rare-variant variance CONDITIONAL on those SNPs.
// (C'C)^+ is formed once per chunk and shared by the target and both flanks.

static void apply_alpha(GenoMat& X, const std::vector<float>& maf, double alpha,
                        const std::vector<float>* wts, int i0) {
    bool unit = (alpha == -1.0);
    double e = (1.0 + alpha) / 2.0;
    int nw = (wts != 0) ? (int) wts->size() : 0;
    for (int j = 0; j < X.cols(); ++j) {
        double w = 1.0;
        if (!unit) { double f = maf[j]; w = (f > 0.0 && f < 1.0) ? std::pow(2.0 * f * (1.0 - f), e) : 1.0; }
        if (nw > 0) { int g = i0 + j; double lw = (g >= 0 && g < nw) ? (double)(*wts)[g] : 1.0;
                      if (lw < 0.0) lw = 0.0; w *= std::sqrt(lw); }
        if (w != 1.0) X.col(j) *= (float) w;
    }
}

// As apply_alpha, but for a non-contiguous SNP selection: the per-SNP weight
// lookup uses sel[j] instead of i0 + j.

// Pascal's triangle up to order MAXCUM, built once. A function-local static has
// thread-safe initialisation in C++11, so this is safe inside the workers.
enum { MAXCUM = 12 };
struct BinomTable {
    double B[MAXCUM + 1][MAXCUM + 1];
    BinomTable() {
        for (int i = 0; i <= MAXCUM; ++i) {
            for (int j = 0; j <= MAXCUM; ++j) B[i][j] = 0.0;
            B[i][0] = 1.0;
            for (int j = 1; j <= i; ++j)
                B[i][j] = B[i-1][j-1] + ((j <= i-1) ? B[i-1][j] : 0.0);
        }
    }
};
static const BinomTable& binoms() { static BinomTable t; return t; }

// central moments -> cumulants (mu[1] must be 0)
static void moms_to_cums(const double* mu, double* kap, int up) {
    const BinomTable& bt = binoms();
    kap[1] = mu[1];
    for (int nn = 2; nn <= up; ++nn) {
        double s = mu[nn];
        for (int mm = 1; mm < nn; ++mm) s -= bt.B[nn-1][mm-1] * kap[mm] * mu[nn-mm];
        kap[nn] = s;
    }
}
// cumulants -> central moments
static void cums_to_moms(const double* kap, double* mu, int up) {
    const BinomTable& bt = binoms();
    mu[0] = 1.0; mu[1] = kap[1];
    for (int nn = 2; nn <= up; ++nn) {
        double s = kap[nn];
        for (int mm = 1; mm < nn; ++mm) s += bt.B[nn-1][mm-1] * kap[mm] * mu[nn-mm];
        mu[nn] = s;
    }
}

// Match kappa1..kappa4 with a two-component spectrum {(l1,v1),(l2,v2)}:
//     kappa_m = 2^{m-1} (m-1)! * sum_j v_j l_j^m
// Damped Newton from several starts. Returns false if no valid (v_j > 0) fit.
// Given the two eigenvalues, the two MULTIPLICITIES solve a 2x2 LINEAR system
// (from kappa1, kappa2), so only a 2-D root find in (l1, l2) remains. Solving it
// that way rather than as a blind 4-D Newton takes the success rate from ~34% to
// essentially 100%, which matters because a failure here silently falls back to
// the weaker variance-only correction.
static bool solve_nu(double a, double b, const double* kq, double* p, double* q) {
    const double det = a * b * (b - a);              // det[[a,b],[a^2,b^2]]
    if (!(std::abs(det) > 1e-300) || !std::isfinite(det)) return false;
    const double r1 = kq[1], r2 = kq[2] / 2.0;
    *p = ( b * b * r1 - b * r2) / det;
    *q = (-a * a * r1 + a * r2) / det;
    return std::isfinite(*p) && std::isfinite(*q);
}
static bool fit_two_component(const double* kq, double* l1, double* v1,
                              double* l2, double* v2) {
    const double sc = std::sqrt(std::abs(kq[2]) / 2.0);
    if (!(sc > 0.0) || !std::isfinite(sc)) return false;
    // residuals in kappa3 / kappa4 as a function of (a,b) alone
    struct R { static bool f(double a, double b, const double* kq,
                             double& r3, double& r4, double& p, double& q) {
        if (!solve_nu(a, b, kq, &p, &q)) return false;
        r3 =  8.0 * (p*a*a*a       + q*b*b*b)       - kq[3];
        r4 = 48.0 * (p*a*a*a*a     + q*b*b*b*b)     - kq[4];
        return std::isfinite(r3) && std::isfinite(r4);
    } };
    // Dense grid of starts. Each solve is a handful of 2x2 operations, so this is
    // free, and it takes the success rate from ~70% to ~100%. Multiple roots
    // exist (the labelling of the two components is arbitrary); any valid one
    // with positive multiplicities reproduces the same distribution.
    static const double GA[6] = {0.25, 0.5, 1.0, 2.0, 4.0, 8.0};
    static const double GB[6] = {-0.25, -0.5, -1.0, -2.0, -4.0, -8.0};
    double best = -1.0;
    for (int st = 0; st < 36; ++st) {
        double a = GA[st / 6]*sc, b = GB[st % 6]*sc;
        for (int it = 0; it < 300; ++it) {
            double r3, r4, p, q;
            if (!R::f(a, b, kq, r3, r4, p, q)) break;
            // Scale the residuals by a NATURAL magnitude, not by kappa_m itself:
            // kappa3 can pass through zero (the two components cancel), and a
            // relative test against it then never converges.
            const double s3 =  8.0 * sc*sc*sc, s4 = 48.0 * sc*sc*sc*sc;
            const double nrm = std::max(std::abs(r3) / std::max(std::abs(kq[3]), s3),
                                        std::abs(r4) / std::max(std::abs(kq[4]), s4));
            if (nrm < 1e-9) {
                if (p > 1e-12 && q > 1e-12 && (best < 0 || nrm < best)) {
                    best = nrm; *l1 = a; *v1 = p; *l2 = b; *v2 = q;
                }
                break;
            }
            // numerical 2x2 Jacobian
            const double ha = 1e-6*std::max(std::abs(a),sc), hb = 1e-6*std::max(std::abs(b),sc);
            double r3a,r4a,r3b,r4b,pp,qq;
            if (!R::f(a+ha,b,kq,r3a,r4a,pp,qq)) break;
            if (!R::f(a,b+hb,kq,r3b,r4b,pp,qq)) break;
            Eigen::Matrix2d J; Eigen::Vector2d F;
            J(0,0)=(r3a-r3)/ha; J(0,1)=(r3b-r3)/hb;
            J(1,0)=(r4a-r4)/ha; J(1,1)=(r4b-r4)/hb;
            F << r3, r4;
            Eigen::Vector2d d = J.fullPivLu().solve(F);
            if (!d.allFinite()) break;
            double damp = 1.0, mx = d.cwiseAbs().maxCoeff();
            if (mx > 0.5*sc) damp = 0.5*sc/mx;          // trust region
            a -= damp*d[0]; b -= damp*d[1];
            if (!std::isfinite(a) || !std::isfinite(b)) break;
            if (std::abs(a-b) < 1e-12*sc) break;
        }
    }
    return best >= 0.0;
}

// ---------------------------------------------------------------------------
// THREE-component fit, matching SIX cumulants.
//
// Truncating at four cumulants is the dominant remaining error: against a Monte
// Carlo reference the four-cumulant fit overshoots the upper tail by ~8% at
// 1e-2 and ~24% at 1e-3, and moving to six roughly halves both. (Replacing the
// saddlepoint with exact characteristic-function inversion, by contrast, changes
// the answer by ~3% -- so the saddlepoint itself is not worth attacking.)
//
// Same structural trick as the two-component case, one dimension up: GIVEN the
// three eigenvalues, the three multiplicities solve a 3x3 LINEAR system taken
// from kappa1..kappa3, leaving only a 3-D root-find on the kappa4..kappa6
// residuals. The system matrix has rows (l_j^m) for m = 1,2,3, whose determinant
// is l1 l2 l3 (l2-l1)(l3-l1)(l3-l2): invertible whenever the values are distinct
// and nonzero.
// ---------------------------------------------------------------------------
static bool solve_nu3(const double* l, const double* kq, double* nu) {
    Eigen::Matrix3d A; Eigen::Vector3d r;
    for (int m = 1; m <= 3; ++m) {
        for (int j = 0; j < 3; ++j) A(m-1, j) = std::pow(l[j], m);
        r[m-1] = kq[m] / (m == 1 ? 1.0 : (m == 2 ? 2.0 : 8.0));
    }
    Eigen::FullPivLU<Eigen::Matrix3d> lu(A);
    if (!lu.isInvertible()) return false;
    Eigen::Vector3d x = lu.solve(r);
    if (!x.allFinite()) return false;
    for (int j = 0; j < 3; ++j) nu[j] = x[j];
    return true;
}
static bool fit_three_component(const double* kq, double* lout, double* nout) {
    const double sc = std::sqrt(std::abs(kq[2]) / 2.0);
    if (!(sc > 0.0) || !std::isfinite(sc)) return false;
    const double w[7] = {0, 1, 2, 8, 48, 384, 3840};      // 2^{m-1}(m-1)!
    struct R { static bool res(const double* l, const double* kq,
                               double* rr, double* nu) {
        if (!solve_nu3(l, kq, nu)) return false;
        const double w[7] = {0, 1, 2, 8, 48, 384, 3840};
        for (int m = 4; m <= 6; ++m) {
            double s = 0.0;
            for (int j = 0; j < 3; ++j) s += nu[j] * std::pow(l[j], m);
            rr[m-4] = w[m] * s - kq[m];
            if (!std::isfinite(rr[m-4])) return false;
        }
        return true;
    } };
    // Starting points. The plain grid alone converges on only about a third of
    // cases; seeding from the (much more reliable) two-component solution, which
    // already matches kappa1..kappa4, roughly doubles that. Those seeds are tried
    // first.
    double S[20][3]; int nS = 0;
    double s1, w1, s2, w2;
    if (fit_two_component(kq, &s1, &w1, &s2, &w2)) {
        static const double thirds[6] = {0.25, 0.5, 2.0, 4.0, -0.5, -2.0};
        for (int q = 0; q < 6; ++q) {
            const double l3 = (thirds[q] > 0 ? s1 : s2) * std::abs(thirds[q]);
            S[nS][0] = s1; S[nS][1] = s2; S[nS][2] = l3; ++nS;
        }
    }
    static const double G[8][3] = {{ 1,-1, 2},{ 2,-1, 0.5},{ 1,-2, 3},
                                   { 3,-0.5,1},{0.5,-3, 2},{ 1,-1, 4},
                                   { 4,-2, 1},{ 2,-3, 0.5}};
    for (int q = 0; q < 8 && nS < 20; ++q) {
        S[nS][0] = G[q][0]*sc; S[nS][1] = G[q][1]*sc; S[nS][2] = G[q][2]*sc; ++nS;
    }
    double best = -1.0;
    for (int st = 0; st < nS; ++st) {
        double l[3] = {S[st][0], S[st][1], S[st][2]};
        for (int it = 0; it < 300; ++it) {
            double rr[3], nu[3];
            if (!R::res(l, kq, rr, nu)) break;
            double nrm = 0.0;
            for (int m = 4; m <= 6; ++m) {
                const double scl = std::max(std::abs(kq[m]), w[m]*std::pow(sc, m));
                nrm = std::max(nrm, std::abs(rr[m-4]) / scl);
            }
            if (nrm < 1e-9) {
                if (nu[0] > 1e-12 && nu[1] > 1e-12 && nu[2] > 1e-12 &&
                    (best < 0 || nrm < best)) {
                    best = nrm;
                    for (int j = 0; j < 3; ++j) { lout[j] = l[j]; nout[j] = nu[j]; }
                }
                break;
            }
            Eigen::Matrix3d J; Eigen::Vector3d F(rr[0], rr[1], rr[2]);
            bool ok = true;
            for (int c = 0; c < 3 && ok; ++c) {
                double lp[3] = {l[0], l[1], l[2]};
                const double h = 1e-6 * std::max(std::abs(l[c]), sc);
                lp[c] += h;
                double r2[3], n2[3];
                if (!R::res(lp, kq, r2, n2)) { ok = false; break; }
                for (int m = 0; m < 3; ++m) J(m, c) = (r2[m] - rr[m]) / h;
            }
            if (!ok) break;
            Eigen::Vector3d d = J.fullPivLu().solve(F);
            if (!d.allFinite()) break;
            double damp = 1.0, mx = d.cwiseAbs().maxCoeff();
            if (mx > 0.5 * sc) damp = 0.5 * sc / mx;
            for (int j = 0; j < 3; ++j) l[j] -= damp * d[j];
            if (!std::isfinite(l[0]) || !std::isfinite(l[1]) || !std::isfinite(l[2])) break;
            if (std::abs(l[0]-l[1]) < 1e-12*sc || std::abs(l[0]-l[2]) < 1e-12*sc ||
                std::abs(l[1]-l[2]) < 1e-12*sc) break;
        }
    }
    return best >= 0.0;
}

struct QuadSpaResult { double p; bool converged; };

// `n_rep` and `nu2` are doubles: the binary 4-cumulant path needs FRACTIONAL
// multiplicities (a scaled chi^2_nu is well defined for any nu > 0). The second
// repeated pair (lam2, nu2) defaults to zero, so the Gaussian path is unchanged.
static void quad_cgf(double t, const std::vector<double>& eig_explicit,
                     double eig_rep, double n_rep, double& K, double& K1, double& K2,
                     double lam2 = 0.0, double nu2 = 0.0,
                     double lam3 = 0.0, double nu3 = 0.0) {
    K = 0.0; K1 = 0.0; K2 = 0.0;
    for (size_t j = 0; j < eig_explicit.size(); ++j) {
        double lam = eig_explicit[j];
        double d = 1.0 - 2.0 * lam * t;
        K  += -0.5 * std::log(d);
        K1 += lam / d;
        K2 += 2.0 * lam * lam / (d * d);
    }
    if (n_rep > 0) {
        double d = 1.0 - 2.0 * eig_rep * t;
        double nr = n_rep;
        K  += -0.5 * nr * std::log(d);
        K1 += nr * eig_rep / d;
        K2 += nr * 2.0 * eig_rep * eig_rep / (d * d);
    }
    if (nu2 > 0) {
        double d = 1.0 - 2.0 * lam2 * t;
        K  += -0.5 * nu2 * std::log(d);
        K1 += nu2 * lam2 / d;
        K2 += nu2 * 2.0 * lam2 * lam2 / (d * d);
    }
    if (nu3 > 0) {
        double d = 1.0 - 2.0 * lam3 * t;
        K  += -0.5 * nu3 * std::log(d);
        K1 += nu3 * lam3 / d;
        K2 += nu3 * 2.0 * lam3 * lam3 / (d * d);
    }
}

static QuadSpaResult quad_spa_solve(
    double s_obs_in, const std::vector<double>& eig_in, double eig_rep_in, double n_rep,
    int max_iter = 100, double tol = 1e-8,
    double lam2_in = 0.0, double nu2 = 0.0,
    double lam3_in = 0.0, double nu3 = 0.0
) {
    QuadSpaResult res; res.p = NA_REAL; res.converged = false;

    // -----------------------------------------------------------------------
    // SCALE NORMALIZATION -- do not remove.
    //
    // The saddlepoint p-value is INVARIANT under a common positive rescaling
    // of (lambda, s_obs): scaling both by k scales the statistic by k and
    // leaves P(X >= s_obs) unchanged (K(t) simply becomes K(kt)). We exploit
    // that to work internally in units where Var = 1.
    //
    // This matters enormously in practice. The eigenvalues here carry the
    // units of a per-window variance component, whose magnitude depends on n
    // and on the phenotype scaling: at biobank n a single 1 Mb window's
    // sigma_hat has SE ~ 1e-7 or smaller, so var0 ~ 1e-14 and the individual
    // lambda are ~1e-9. Every hard-coded tolerance below (convergence
    // threshold, degeneracy floor, domain epsilon, the |t| ~ 0 test) would
    // then be comparing against constants many orders of magnitude too large,
    // and the solver would bail out on EVERY window -- returning NA
    // everywhere, while working fine on smaller test data where the same
    // quantities happen to be O(1e-3). That was a real bug.
    //
    // After normalizing, var0 == 1 exactly, mean0 ~ 0, s_obs is a z-score,
    // |lambda| <= 1/sqrt(2), and t is O(1) -- so the dimensionless constants
    // below are meaningful regardless of n or phenotype scaling.
    // -----------------------------------------------------------------------
    double mean_raw = 0.0, var_raw = 0.0;
    for (size_t j = 0; j < eig_in.size(); ++j) {
        double lam = eig_in[j];
        mean_raw += lam; var_raw += 2.0 * lam * lam;
    }
    if (n_rep > 0) {
        mean_raw += n_rep * eig_rep_in;
        var_raw  += n_rep * 2.0 * eig_rep_in * eig_rep_in;
    }
    if (nu2 > 0) {
        mean_raw += nu2 * lam2_in;
        var_raw  += nu2 * 2.0 * lam2_in * lam2_in;
    }
    if (nu3 > 0) {
        mean_raw += nu3 * lam3_in;
        var_raw  += nu3 * 2.0 * lam3_in * lam3_in;
    }
    // Only a genuinely degenerate (zero / non-finite) variance is unusable.
    if (!(var_raw > 0.0) || !std::isfinite(var_raw)) return res;

    const double sd_raw = std::sqrt(var_raw);
    const double scale  = 1.0 / sd_raw;

    std::vector<double> eig_explicit(eig_in.size());
    for (size_t j = 0; j < eig_in.size(); ++j) eig_explicit[j] = eig_in[j] * scale;
    const double eig_rep = eig_rep_in * scale;
    const double lam2    = lam2_in    * scale;
    const double lam3    = lam3_in    * scale;
    const double s_obs   = s_obs_in   * scale;

    double lam_min = (n_rep > 0) ? eig_rep : 0.0, lam_max = lam_min;
    for (size_t j = 0; j < eig_explicit.size(); ++j) {
        double lam = eig_explicit[j];
        if (lam < lam_min) lam_min = lam;
        if (lam > lam_max) lam_max = lam;
    }
    if (nu2 > 0) { if (lam2 < lam_min) lam_min = lam2; if (lam2 > lam_max) lam_max = lam2; }
    if (nu3 > 0) { if (lam3 < lam_min) lam_min = lam3; if (lam3 > lam_max) lam_max = lam3; }
    // Dimensionless now: |lambda| <= 1/sqrt(2), so 1e-14 simply means
    // "no positive (negative) eigenvalue, hence no bound on t in that direction".
    const double lam_eps = 1e-14;
    double t_hi = (lam_max >  lam_eps) ? (1.0 / (2.0 * lam_max)) :  1e12;
    double t_lo = (lam_min < -lam_eps) ? (1.0 / (2.0 * lam_min)) : -1e12;
    double margin = 1e-6;
    double t_hi_safe = (lam_max >  lam_eps) ? t_hi * (1.0 - margin) : t_hi;
    double t_lo_safe = (lam_min < -lam_eps) ? t_lo * (1.0 - margin) : t_lo;

    double t = 0.0, K, K1, K2;
    quad_cgf(0.0, eig_explicit, eig_rep, n_rep, K, K1, K2, lam2, nu2, lam3, nu3);
    double mean0 = K1, var0 = K2;          // var0 == 1 up to rounding
    if (!(var0 > 0.0)) return res;

    bool converged = false;
    for (int it = 0; it < max_iter; ++it) {
        quad_cgf(t, eig_explicit, eig_rep, n_rep, K, K1, K2, lam2, nu2, lam3, nu3);
        double diff = K1 - s_obs;
        if (std::abs(diff) < tol * std::max(1.0, std::abs(s_obs - mean0))) { converged = true; break; }
        if (K2 < 1e-14) break;
        double step = diff / K2;
        double t_new = t - step;
        int halvings = 0;
        while ((t_new <= t_lo_safe || t_new >= t_hi_safe) && halvings < 40) {
            step *= 0.5; t_new = t - step; ++halvings;
        }
        if (t_new <= t_lo_safe || t_new >= t_hi_safe) break;   // safeguard exhausted
        t = t_new;
    }
    if (!converged) return res;

    quad_cgf(t, eig_explicit, eig_rep, n_rep, K, K1, K2, lam2, nu2, lam3, nu3);
    if (K2 <= 0) return res;

    // Removable singularity at t ~ 0: the normal approximation is exact in
    // that limit, so use it directly rather than divide by ~0.
    if (std::abs(t) < 1e-7) {
        double z = (s_obs - mean0) / std::sqrt(var0);
        res.p = std::erfc(std::abs(z) / std::sqrt(2.0));
        res.converged = true;
        return res;
    }

    double w = ((t > 0) ? 1.0 : -1.0) * std::sqrt(std::max(0.0, 2.0 * (t * s_obs - K)));
    double u = t * std::sqrt(K2);
    double Phi_w = 0.5 * std::erfc(-w / std::sqrt(2.0));
    double phi_w = std::exp(-0.5 * w * w) / std::sqrt(2.0 * M_PI);

    // Upper tail computed DIRECTLY as 0.5*erfc(w/sqrt2). Writing it as
    // 1 - Phi(w) cancels catastrophically: past w ~ 8.3 it returns exactly 0,
    // which is where p_spa = 0 comes from. The lower branch has no cancellation.
    double p_one = (w >= 0)
        ? 0.5 * std::erfc(w / std::sqrt(2.0)) + phi_w * (1.0 / u - 1.0 / w)
        : Phi_w - phi_w * (1.0 / u - 1.0 / w);
    if (p_one < 0.0) p_one = 0.0;
    if (p_one > 1.0) p_one = 1.0;

    res.p = std::min(1.0, 2.0 * p_one);
    res.converged = true;
    return res;
}





// ===========================================================================
// BINARY CO-HERITABILITY:  exact conditional saddlepoint
// ===========================================================================
// Q = y_1' M_c y_2 is LINEAR in y_2. Hold trait 1 fixed, put w = M_c y_1, and
//     Q = w' y_2
// is a weighted sum of the second trait's 0/1 outcomes -- which has an EXACT
// closed-form cumulant generating function. No chi-square mixture, no cumulant
// truncation, no normality. This is the Dey/SAIGE saddlepoint, and it is the
// reason the bivariate binary problem is EASIER than the univariate one:
// y' M y has the same y on both sides and cannot be linearised.
//
// The analysis phenotype is the covariate residual of the 0/1 trait, rescaled:
//     y_i = (b_i - yhat_i) / s ,   b_i in {0,1}
// so it still takes exactly TWO values per individual, a gap of delta = 1/s
// apart. Writing y_i = ylo_i + delta * b_i,
//     Q = A + sum_i c_i b_i ,   A = w' ylo ,  c_i = delta * w_i
// with b_i ~ Bernoulli(pi_i) independent given trait 1. Hence
//     K(t) = A t + sum_i log(1 - pi_i + pi_i e^{c_i t})
// exactly.
//
// s is recovered without extra input: y = (I-H)b/s and b'(I-H)b = ||(I-H)b||^2,
// so b'y = (n-1)s, i.e. delta = (n-1) / (b'y).
//
// Validated by simulation against the Gaussian spectrum SPA. At two binary
// traits with prevalence 2% and correlation 0.3 the Gaussian null gives
// 0.0994 / 0.0384 / 0.0156 at alpha = 0.05 / 0.01 / 0.001 (a 16-fold excess in
// the tail); the conditional Bernoulli saddlepoint gives 0.0498 / 0.0096 /
// 0.0012.
// ===========================================================================
static QuadSpaResult bern_spa_solve(double q_obs, const std::vector<double>& cc,
                                    const std::vector<float>& pi, double A,
                                    int max_iter = 100, double tol = 1e-10) {
    QuadSpaResult res; res.p = NA_REAL; res.converged = false;
    const int n = (int) cc.size();
    double mean0 = A, var0 = 0.0;
    for (int i = 0; i < n; ++i) {
        const double pv = pi[i];
        mean0 += cc[i] * pv;
        var0  += cc[i] * cc[i] * pv * (1.0 - pv);
    }
    if (!(var0 > 0.0) || !std::isfinite(var0)) return res;
    // Same scale normalisation as quad_spa_solve: work in units where Var = 1,
    // otherwise every tolerance below is compared against the wrong magnitude.
    const double sd = std::sqrt(var0), sc = 1.0 / sd;
    std::vector<double> c(n);
    for (int i = 0; i < n; ++i) c[i] = cc[i] * sc;
    const double An = A * sc, qn = q_obs * sc;

    double t = 0.0, Kv = 0.0, K1 = 0.0, K2 = 0.0;
    bool conv = false;
    for (int it = 0; it < max_iter; ++it) {
        Kv = An * t; K1 = An; K2 = 0.0;
        for (int i = 0; i < n; ++i) {
            const double z = c[i] * t, pv = pi[i];
            const double m = (z > 0.0) ? z : 0.0;          // log-sum-exp shift
            const double e0 = (1.0 - pv) * std::exp(-m), e1 = pv * std::exp(z - m);
            const double den = e0 + e1;
            Kv += m + std::log(den);
            const double r = e1 / den;
            K1 += c[i] * r;
            K2 += c[i] * c[i] * r * (1.0 - r);
        }
        const double diff = K1 - qn;
        if (std::abs(diff) < tol * std::max(1.0, std::abs(qn - An))) { conv = true; break; }
        if (!(K2 > 1e-14)) break;
        double step = diff / K2;
        if (step >  2.0) step =  2.0;                       // the domain is all of R,
        if (step < -2.0) step = -2.0;                       // but keep exp() in range
        t -= step;
    }
    if (!conv) return res;
    if (!(K2 > 0.0)) return res;

    if (std::abs(t) < 1e-9) {                               // removable singularity
        const double z = (qn - (An + K1 - An)) / std::sqrt(K2);
        res.p = std::erfc(std::abs(z) / std::sqrt(2.0));
        res.converged = true; return res;
    }
    const double w = ((t > 0) ? 1.0 : -1.0)
                   * std::sqrt(std::max(0.0, 2.0 * (t * qn - Kv)));
    const double u = t * std::sqrt(K2);
    const double Phi_w = 0.5 * std::erfc(-w / std::sqrt(2.0));
    const double phi_w = std::exp(-0.5 * w * w) / std::sqrt(2.0 * M_PI);
    // see quad_spa_solve: the upper tail must not be written as 1 - Phi(w)
    double p_one = (w >= 0) ? 0.5 * std::erfc(w / std::sqrt(2.0))
                              + phi_w * (1.0 / u - 1.0 / w)
                            : Phi_w - phi_w * (1.0 / u - 1.0 / w);
    if (p_one < 0.0) p_one = 0.0;
    if (p_one > 1.0) p_one = 1.0;
    res.p = std::min(1.0, 2.0 * p_one);
    res.converged = true;
    return res;
}

static void coher_mdiag(const GenoMat& V, const std::vector<int>& off, int C,
                        const std::vector<double>& g, const Eigen::MatrixXd& Tinv,
                        int c, int env, int n, Eigen::VectorXd& md) {
    md = Eigen::VectorXd::Constant(n, Tinv(c, env));
    for (int a = 0; a < C; ++a) {
        const double dm = Tinv(c, a) * g[a];
        if (dm == 0.0) continue;
        for (int j = off[a]; j < off[a + 1]; ++j)
            md.array() += dm * V.col(j).cast<double>().array().square();
    }
}

// w = M_c y, from the chunk's V and the already-computed u = V' y. One GEMV.
static void coher_wvec(const GenoMat& V, const std::vector<int>& off, int C, int K,
                       const std::vector<double>& g, const Eigen::MatrixXd& Tinv,
                       int c, int env, const Eigen::VectorXd& u,
                       const Eigen::VectorXd& y, Eigen::VectorXd& w) {
    GenoMat d(K, 1);
    for (int a = 0; a < C; ++a) {
        const double dm = Tinv(c, a) * g[a];
        for (int j = off[a]; j < off[a + 1]; ++j) d(j, 0) = (float)(dm * u[j]);
    }
    w = (V * d).col(0).cast<double>() + Tinv(c, env) * y;
}

// Per-pair inputs for the exact conditional saddlepoint (binary coher).
struct CoherBin {
    std::vector<double> k22;                 // joint 4th cumulant of the pair
    std::vector<int>    cond_t, rand_t;      // fixed trait / random trait
    std::vector<double> cdelta;              // gap between the two phenotype values
    std::vector< std::vector<float> > cpi;   // P(random = case | fixed trait)
    std::vector< std::vector<float> > cylo;  // phenotype value if a control
    bool active;
    CoherBin() : active(false) {}
};

// ===========================================================================
// CO-HERITABILITY  (coher = TRUE)
// ===========================================================================
// The estimator is the SAME linear system. Because sigma_hat = T^+ q is linear
// in q, feeding it the CROSS moments
//     q_a = g_a * sum_{j in a} (V' y1)_j (V' y2)_j ,   q_env = y1' y2
// returns the genetic COVARIANCE of each component instead of the variance.
// T is untouched: it depends only on the kernels, not on the phenotype.
//
// The null is where the two differ. Writing Q = y1' M_c y2 and stacking
// z = (y1, y2), Q = z' B z with B = [[0, M_c/2],[M_c/2, 0]] and
// Omega = [[S1, S12],[S12, S2]]. B is INDEFINITE, so the null is a chi-square
// mixture with eigenvalues of BOTH signs -- exactly the regime where normal and
// Satterthwaite approximations fail in both tails and where the saddlepoint
// earns its keep. quad_spa_solve already brackets t on two sides and already
// returns a two-sided p-value, so it is reused unchanged.
//
// H0 IS NOT THE UNIVARIATE H0. Under "no genetic covariance at this locus" the
// two traits may each still be heritable here, so the tested component is KEPT
// in S1 and S2 and nulled only in S12. (The univariate test nulls it in Sigma_0
// itself.) Getting this wrong makes the test conservative at heritable loci.
//
// Exact variance, from tr((B Omega)^2):
//     Var(Q) = tr(M S1 M S2) + tr((M S12)^2)
// which in the K-space representation used throughout is
//     tr(A S1 A S2) + tr((A S12)^2) + (n-K) c_env^2 (e1 e2 + e12^2).
//
// On col(V)^perp all four matrices act as scalars, so those n-K directions each
// contribute a 2x2 problem with eigenvalues (c_env/2)(e12 +/- sqrt(e1 e2)) --
// one positive, one negative. They go into quad_spa_solve's repeated slots.
//
// Validated against a direct 2n-dimensional simulation: the reduced spectrum
// matches the full one to 7 digits, the variance formula matches the empirical
// variance to 0.4%, and the SPA p-values are uniform (0.049 / 0.0095 / 0.0014
// observed at alpha = 0.05 / 0.01 / 0.001).
// ===========================================================================
struct CoherFit {
    Eigen::MatrixXd S1, S2, S12, Ac;
    double var0, e1, e2, e12, c_env;
    bool ok;
    CoherFit() : var0(0), e1(0), e2(0), e12(0), c_env(0), ok(false) {}
};

// S = L' D L + e I for a block-constant D given by per-component values.
static void coher_block(const Eigen::MatrixXd& L,
                        const std::vector<Eigen::MatrixXd>& W, bool use_wcache,
                        const std::vector<int>& off, int C, int K,
                        const std::vector<double>& d, double e,
                        Eigen::MatrixXd& out) {
    if (use_wcache) {
        out.setZero(K, K);
        for (int a = 0; a < C; ++a) if (d[a] != 0.0) out.noalias() += d[a] * W[a];
    } else {
        Eigen::VectorXd D(K);
        for (int a = 0; a < C; ++a)
            for (int j = off[a]; j < off[a + 1]; ++j) D[j] = d[a];
        out.noalias() = L.transpose() * (D.asDiagonal() * L);
    }
    out.diagonal().array() += e;
}

// Null variance of Q = y1' M_c y2 under H0: cov(component c) = 0.
static void coher_fit(int c, const Eigen::MatrixXd& Tinv, int C, int env,
                      const std::vector<double>& g, const std::vector<int>& off,
                      int K, int n, const Eigen::MatrixXd& L,
                      const std::vector<Eigen::MatrixXd>& W, bool use_wcache,
                      const Eigen::MatrixXd& Ac,
                      const Eigen::VectorXd& sg1, const Eigen::VectorXd& sg2,
                      const Eigen::VectorXd& s12, double Vp1, double Vp2,
                      double kappa22, double Sb2c, CoherFit& F) {
    std::vector<double> d1(C), d2(C), d12(C);
    for (int a = 0; a < C; ++a) {
        // DO NOT CLAMP THE TESTED COMPONENT. It is a nuisance here, not the
        // null hypothesis, so it must stay in Sigma_1 and Sigma_2 -- zeroing it
        // gives 84% type I error at loci where the traits really are heritable.
        // But clamping it at zero is just as wrong in the other direction: where
        // the true target h2 is 0, sigma_hat is symmetric about 0, the clamp
        // fires ~half the time and returns E[max(0,sigma_hat)] ~ 0.4 sd instead
        // of 0 -- for BOTH traits, multiplied together inside
        // tr(M Sigma_1 M Sigma_2). Measured effect: SE inflated ~15%, variance
        // ~33%, lambda_GC ~ 0.8, i.e. uniformly DEFLATED p-values with a
        // straight QQ line of slope ~0.8.
        // The clamp is kept on the background components, whose true values sit
        // well away from zero so it essentially never fires, and where it does
        // useful work keeping Omega positive definite.
        d1[a]  = (a == c ? sg1[a] : (sg1[a] > 0.0 ? sg1[a] : 0.0)) * g[a];
        d2[a]  = (a == c ? sg2[a] : (sg2[a] > 0.0 ? sg2[a] : 0.0)) * g[a];
        d12[a] = s12[a] * g[a];
    }
    d12[c] = 0.0;                                        // the tested co-component
    F.e1  = std::max(sg1[env], 1e-8 * (Vp1 > 0 ? Vp1 : 1.0));
    F.e2  = std::max(sg2[env], 1e-8 * (Vp2 > 0 ? Vp2 : 1.0));
    F.e12 = s12[env];
    F.c_env = Tinv(c, env);
    // |e12| must respect Cauchy-Schwarz or Omega is not a covariance matrix.
    const double ecap = 0.999 * std::sqrt(F.e1 * F.e2);
    if (F.e12 >  ecap) F.e12 =  ecap;
    if (F.e12 < -ecap) F.e12 = -ecap;

    coher_block(L, W, use_wcache, off, C, K, d1,  F.e1,  F.S1);
    coher_block(L, W, use_wcache, off, C, K, d2,  F.e2,  F.S2);
    coher_block(L, W, use_wcache, off, C, K, d12, F.e12, F.S12);
    F.Ac = Ac;

    Eigen::MatrixXd N1, N2, N12;
    N1.noalias()  = F.Ac * F.S1;
    N2.noalias()  = F.Ac * F.S2;
    N12.noalias() = F.Ac * F.S12;
    double v = (N1.cwiseProduct(N2.transpose())).sum()
             + (N12.cwiseProduct(N12.transpose())).sum();
    const long n_rep = (long) n - K;
    if (n_rep > 0)
        v += (double) n_rep * F.c_env * F.c_env * (F.e1 * F.e2 + F.e12 * F.e12);
    // BINARY CORRECTION. Working with M directly (not the stacked B, whose
    // diagonal is zero) the exact variance is
    //     Var(y1' M y2) = (v1 v2 + v12^2) tr(M^2) + kappa_22 * sum_i M_ii^2
    // which collapses to the univariate 2 v^2 tr(M^2) + kappa_4 sum M_ii^2 when
    // y1 = y2. So the diagonal term survives, with the JOINT fourth cumulant in
    // place of kappa_4, and it reuses the same sum_i M_ii^2 (Sb2).
    // Note kappa_22 = 0 when the traits are independent, whatever their
    // marginals: binarity enters a covariance null only through trait-trait
    // dependence, and only appreciably when BOTH traits are rare.
    if (kappa22 != 0.0 && Sb2c > 0.0) v += kappa22 * Sb2c;
    F.var0 = v;
    F.ok = (v > 0.0) && std::isfinite(v);
}

// Full null spectrum: 2K explicit eigenvalues plus two repeated groups.
static bool coher_spectrum(const CoherFit& F, int K, int n,
                           std::vector<double>& eig,
                           double& eig_rep, double& lam2, double& n_rep) {
    const int K2 = 2 * K;
    Eigen::MatrixXd Om(K2, K2);
    Om.topLeftCorner(K, K)     = F.S1;
    Om.topRightCorner(K, K)    = F.S12;
    Om.bottomLeftCorner(K, K)  = F.S12;
    Om.bottomRightCorner(K, K) = F.S2;
    // Plug-in estimates need not give a PSD Omega. Shrink the cross block until
    // they do; shrinking S12 toward 0 moves toward independence, the
    // conservative direction for a covariance test.
    // An unclamped negative estimate for the tested component can, rarely, tip a
    // diagonal block indefinite. Shrinking the cross block cannot fix that, so
    // pull the diagonal blocks toward their environment-only form first.
    Eigen::LLT<Eigen::MatrixXd> llt;
    double shrink = 1.0, dshrink = 1.0; bool have = false;
    for (int it = 0; it < 16; ++it) {
        llt.compute(Om);
        if (llt.info() == Eigen::Success) { have = true; break; }
        if (Eigen::LLT<Eigen::MatrixXd>(Om.topLeftCorner(K, K)).info() != Eigen::Success ||
            Eigen::LLT<Eigen::MatrixXd>(Om.bottomRightCorner(K, K)).info() != Eigen::Success) {
            dshrink *= 0.5;
            Om.topLeftCorner(K, K) =
                dshrink * F.S1 + (1.0 - dshrink) * F.e1 * Eigen::MatrixXd::Identity(K, K);
            Om.bottomRightCorner(K, K) =
                dshrink * F.S2 + (1.0 - dshrink) * F.e2 * Eigen::MatrixXd::Identity(K, K);
        } else {
            shrink *= 0.5;
        }
        Om.topRightCorner(K, K)   = shrink * F.S12;
        Om.bottomLeftCorner(K, K) = shrink * F.S12;
    }
    if (!have) return false;
    Eigen::MatrixXd Cf = llt.matrixL();
    // Bk = [[0, A/2],[A/2, 0]];  Asym = Cf' Bk Cf, symmetric, 2K x 2K.
    Eigen::MatrixXd BC(K2, K2);
    BC.topRows(K).noalias()    = 0.5 * (F.Ac * Cf.bottomRows(K));
    BC.bottomRows(K).noalias() = 0.5 * (F.Ac * Cf.topRows(K));
    Eigen::MatrixXd Asym; Asym.noalias() = Cf.transpose() * BC;
    Asym = (0.5 * (Asym + Asym.transpose())).eval();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ses(Asym, Eigen::EigenvaluesOnly);
    if (ses.info() != Eigen::Success) return false;
    const Eigen::VectorXd& ev = ses.eigenvalues();
    eig.assign(ev.data(), ev.data() + K2);
    const double root = std::sqrt(std::max(0.0, F.e1 * F.e2));
    eig_rep = 0.5 * F.c_env * (F.e12 + root);
    lam2    = 0.5 * F.c_env * (F.e12 - root);
    n_rep   = (double) (n - K);
    return true;
}

// ===========================================================================
// Context
// ===========================================================================
// LDAK-KVIK step-1 LOCO PRS, one source (e.g. the common-SNP panel). The file
// is  FID IID Chr1 ... Chr22 , each column the PRS built from every chromosome
// EXCEPT that one. Because a single fit is partitioned by chromosome,
//     loco_c = full - part_c   =>   sum_c loco_c = (nchr - 1) * full
// so the FULL PRS is recovered exactly as sum_c loco_c / (nchr - 1), and the
// per-chromosome contribution as full - loco_c. Verify this once against
// `--calc-scores` on the .effects file; if KVIK ever refits per fold the
// identity becomes approximate.
struct PrsSource {
    std::string path;
    std::vector<std::string> chr_names;          // as in the file, "Chr1" -> "1"
    std::vector<int> chr_index;                  // column -> ctx.chr_order index, -1
    Eigen::MatrixXd loco;                        // n x nchr, analysis order
    Eigen::VectorXd full;                        // n
};

struct ChunkContext {
    std::string prefix; int n_total, n_snps; BimInfo bim;
    std::vector<int> geno_keep;
    std::vector<int> pheno_keep;                 // analysis row -> phenotype-matrix row
    Eigen::MatrixXd Y;                           // input phenotype, analysis order
    CharacterVector trait_names;
    int n_inds, n_pheno;
    std::vector<std::string> chr_order;
    std::vector<int> chr_lo, chr_hi;
    double alpha;
    GenoMat covZ, covM;
    std::vector<std::string> analysis_iid_s;
    // every SNP belongs to exactly one category; -1 = no category (background)
    int n_cat; std::vector<std::string> cat_names; std::vector<int> snp_cat;
    BedReader bed;
    // PRS sources and the per-category mask used in pass 2
    std::vector<PrsSource> prs;
    std::vector< std::vector<unsigned char> > prs_mask;   // [cat][source]
    // binary raw 0/1 per trait (only when binary && coher), analysis order
    std::vector< std::vector<unsigned char> > braw;
    std::vector<double> prev, bdelta;
};

static std::string strip_chr(const std::string& s) {
    if (s.size() > 3 && (s.compare(0, 3, "Chr") == 0 || s.compare(0, 3, "chr") == 0 ||
                         s.compare(0, 3, "CHR") == 0)) return s.substr(3);
    return s;
}

static void prs_load(const std::string& path, ChunkContext& ctx, PrsSource& S) {
    S.path = path;
    std::ifstream f(path.c_str());
    if (!f.is_open()) stop("Cannot open LOCO PRS file: " + path);
    std::string line;
    if (!std::getline(f, line)) stop("Empty LOCO PRS file: " + path);
    std::istringstream hs(line); std::string tok; std::vector<std::string> hdr;
    while (hs >> tok) hdr.push_back(tok);
    if (hdr.size() < 3) stop("LOCO PRS header needs FID IID Chr...: " + path);
    const int nchr = (int) hdr.size() - 2;
    S.chr_names.resize(nchr); S.chr_index.assign(nchr, -1);
    std::map<std::string,int> chrpos;
    for (size_t i = 0; i < ctx.chr_order.size(); ++i)
        chrpos[strip_chr(ctx.chr_order[i])] = (int) i;
    for (int c = 0; c < nchr; ++c) {
        S.chr_names[c] = strip_chr(hdr[c + 2]);
        std::map<std::string,int>::const_iterator it = chrpos.find(S.chr_names[c]);
        if (it != chrpos.end()) S.chr_index[c] = it->second;
    }
    std::map<std::string,int> want;
    for (int i = 0; i < ctx.n_inds; ++i) want[ctx.analysis_iid_s[i]] = i;
    S.loco.setConstant(ctx.n_inds, nchr, std::numeric_limits<double>::quiet_NaN());
    long nread = 0;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        std::istringstream ls(line); std::string fid, iid;
        if (!(ls >> fid >> iid)) continue;
        std::map<std::string,int>::const_iterator it = want.find(iid);
        if (it == want.end()) it = want.find(fid);
        if (it == want.end()) continue;
        const int r = it->second;
        for (int c = 0; c < nchr; ++c) { double v; if (!(ls >> v)) stop("Short row in " + path); S.loco(r, c) = v; }
        ++nread;
    }
    for (int i = 0; i < ctx.n_inds; ++i)
        if (std::isnan(S.loco(i, 0)))
            stop("LOCO PRS file " + path + " is missing analysis individual " + ctx.analysis_iid_s[i]);
    S.full = S.loco.rowwise().sum() / (double)(nchr - 1);
    int matched = 0; for (int c = 0; c < nchr; ++c) if (S.chr_index[c] >= 0) ++matched;
    Rcout << "LOCO PRS " << path << ": " << nchr << " chromosome columns, " << nread
          << " rows matched, " << matched << " columns map onto the .bim chromosomes.\n"
          << "  full PRS = sum(LOCO)/" << (nchr - 1) << ";  sd(full) per column:";
    // sanity: sd of the recovered full PRS
    double m = S.full.mean(), v = (S.full.array() - m).square().sum() / (ctx.n_inds - 1);
    Rcout << " " << std::sqrt(v) << "\n";
}

// Phenotype variant + everything the tester needs that depends on it.
// Built once per chromosome (pass 1: y - sum LOCO_c) and once per distinct
// prs_mask row (pass 2: y - sum full).
struct PhenoStats {
    Eigen::MatrixXd Y; GenoMat Yf;
    std::vector<double> Vp, yty, kur, ycross;
    CoherBin cb;
};

static ChunkContext setup_context(const std::string& filename, const SEXP pheno_mat,
                                  double alpha,
                                  Rcpp::Nullable<Rcpp::NumericMatrix> covariates,
                                  Rcpp::Nullable<Rcpp::IntegerMatrix> annotation,
                                  Rcpp::Nullable<Rcpp::CharacterVector> annot_names) {
    ChunkContext ctx;
    ctx.prefix = filename;
    ctx.n_snps = count_lines(filename + ".bim");
    List fam = read_fam_file(filename);
    CharacterVector geno_iid = fam["iid"];
    ctx.n_total = geno_iid.size();
    ctx.alpha = alpha;
    Rcout << "Genotype file: " << ctx.n_total << " individuals, " << ctx.n_snps << " SNPs\n";
    ctx.bim = read_bim_positions(filename + ".bim");
    if (ctx.bim.n_snps != ctx.n_snps) stop("Mismatch between .bim line count and parsed positions");
    calibrate_bed_reader(ctx.bed, filename, ctx.n_total, ctx.n_snps, "genotype");

    // categories: one per annotation column; a SNP with an all-zero row is
    // background only. No annotation = one category holding every SNP.
    ctx.snp_cat.assign(ctx.n_snps, -1);
    if (annotation.isNotNull()) {
        Rcpp::IntegerMatrix am(annotation.get());
        if (am.nrow() != ctx.n_snps) stop("annotation must have one row per SNP in the .bim file");
        ctx.n_cat = am.ncol();
        if (ctx.n_cat < 1) stop("annotation needs at least one column");
        for (int j = 0; j < ctx.n_snps; ++j)
            for (int c = 0; c < ctx.n_cat; ++c) if (am(j, c) == 1) { ctx.snp_cat[j] = c; break; }
        if (annot_names.isNotNull()) {
            Rcpp::CharacterVector nm(annot_names.get());
            if (nm.size() != ctx.n_cat) stop("annot_names length must equal the number of annotation columns");
            for (int c = 0; c < ctx.n_cat; ++c) ctx.cat_names.push_back(as<std::string>(nm[c]));
        } else for (int c = 0; c < ctx.n_cat; ++c) ctx.cat_names.push_back("cat" + std::to_string(c + 1));
    } else {
        ctx.n_cat = 1; ctx.cat_names.push_back("ALL");
        for (int j = 0; j < ctx.n_snps; ++j) ctx.snp_cat[j] = 0;
    }
    {
        std::vector<long> cnt(ctx.n_cat, 0); long bg = 0;
        for (int j = 0; j < ctx.n_snps; ++j) { if (ctx.snp_cat[j] >= 0) ++cnt[ctx.snp_cat[j]]; else ++bg; }
        Rcout << "Categories:";
        for (int c = 0; c < ctx.n_cat; ++c) Rcout << " " << ctx.cat_names[c] << "=" << cnt[c];
        if (bg > 0) Rcout << "  (uncategorised, background only: " << bg << ")";
        Rcout << "\n";
    }

    Rcpp::NumericMatrix pheno; CharacterVector pheno_ids;
    if (Rf_isMatrix(pheno_mat) && !Rf_isNull(rownames(pheno_mat))) {
        pheno = as<NumericMatrix>(pheno_mat);
        pheno_ids = rownames(pheno_mat);
        SEXP cn = colnames(pheno_mat);
        if (!Rf_isNull(cn)) ctx.trait_names = cn;
    } else stop("Phenotype must be a numeric matrix with IDs as rownames");
    ctx.n_pheno = pheno.cols();
    if (ctx.trait_names.size() != ctx.n_pheno) {
        ctx.trait_names = CharacterVector(ctx.n_pheno);
        for (int j = 0; j < ctx.n_pheno; ++j) ctx.trait_names[j] = "trait" + std::to_string(j + 1);
    }
    IntegerVector match_idx = match(geno_iid, pheno_ids);
    for (int i = 0; i < match_idx.size(); ++i)
        if (match_idx[i] != NA_INTEGER) { ctx.geno_keep.push_back(i); ctx.pheno_keep.push_back(match_idx[i] - 1); }
    if (ctx.geno_keep.empty()) stop("No overlapping individuals between genotype and phenotype");
    ctx.n_inds = (int) ctx.geno_keep.size();
    Rcout << "Individuals with complete data: " << ctx.n_inds << "\n";
    ctx.Y.resize(ctx.n_inds, ctx.n_pheno);
    for (int i = 0; i < ctx.n_inds; ++i)
        for (int j = 0; j < ctx.n_pheno; ++j) ctx.Y(i, j) = pheno(ctx.pheno_keep[i], j);

    if (covariates.isNotNull()) {
        Rcpp::NumericMatrix cv(covariates.get());
        if (cv.nrow() != pheno.nrow()) stop("covariates must have one row per row of the phenotype matrix");
        int qc = cv.ncol();
        if (qc > 0) {
            Eigen::MatrixXd Zd(ctx.n_inds, qc);
            for (int i = 0; i < ctx.n_inds; ++i)
                for (int j = 0; j < qc; ++j) { double v = cv(ctx.pheno_keep[i], j); Zd(i, j) = ISNAN(v) ? 0.0 : v; }
            for (int j = 0; j < qc; ++j) Zd.col(j).array() -= Zd.col(j).mean();
            Eigen::MatrixXd ZtZ = Zd.transpose() * Zd;
            Eigen::MatrixXd Md = ZtZ.completeOrthogonalDecomposition().pseudoInverse() * Zd.transpose();
            ctx.covZ = Zd.cast<float>(); ctx.covM = Md.cast<float>();
            Rcout << "Regressing " << qc << " covariates out of genotypes\n";
        }
    }
    ctx.analysis_iid_s.resize(ctx.n_inds);
    for (int i = 0; i < ctx.n_inds; ++i)
        ctx.analysis_iid_s[i] = Rcpp::as<std::string>(geno_iid[ctx.geno_keep[i]]);
    for (int j = 0; j < ctx.n_snps; ++j) {
        if (ctx.chr_order.empty() || ctx.bim.chr[j] != ctx.chr_order.back()) {
            ctx.chr_order.push_back(ctx.bim.chr[j]);
            ctx.chr_lo.push_back(j); ctx.chr_hi.push_back(j);
        } else ctx.chr_hi.back() = j;
    }
    return ctx;
}

// Build a phenotype variant: Yv = Y - offset (n x P), plus every y-dependent
// quantity the tester consumes. Vp here is the VARIANT's variance -- used for
// the SPA null floors and the binary cumulant path; the h2 columns are divided
// by the ORIGINAL phenotypic variance in the driver, so heritabilities from
// both passes are on one scale.
static void build_pheno_stats(const ChunkContext& ctx, const Eigen::MatrixXd& offset,
                              bool binary, bool coher,
                              const std::vector< std::pair<int,int> >& pairs,
                              PhenoStats& S) {
    const int n = ctx.n_inds, P = ctx.n_pheno;
    S.Y = ctx.Y;
    if (offset.rows() == n) for (int t = 0; t < P; ++t) S.Y.col(t) -= offset.col(std::min(t, (int) offset.cols() - 1));
    S.Yf = S.Y.cast<float>();
    S.Vp.assign(P, 0.0); S.yty.assign(P, 0.0);
    for (int t = 0; t < P; ++t) {
        double m = S.Y.col(t).mean();
        S.Vp[t]  = (S.Y.col(t).array() - m).square().sum() / (n - 1);
        S.yty[t] = S.Y.col(t).squaredNorm();
    }
    S.kur.clear();
    if (binary) {
        S.kur.assign((size_t) P * (MAXCUM + 1), 0.0);
        for (int t = 0; t < P; ++t) {
            Eigen::ArrayXd d = S.Y.col(t).array() - S.Y.col(t).mean();
            double mu[MAXCUM + 1]; mu[0] = 1.0; mu[1] = 0.0;
            Eigen::ArrayXd pw = d;
            for (int r = 2; r <= MAXCUM; ++r) { pw = pw * d; mu[r] = pw.sum() / (double) n; }
            double kp[MAXCUM + 1]; moms_to_cums(mu, kp, MAXCUM);
            for (int r = 1; r <= MAXCUM; ++r) S.kur[(size_t) t * (MAXCUM + 1) + r] = kp[r];
        }
    }
    S.ycross.clear(); S.cb = CoherBin();
    if (!coher) return;
    const int NP = (int) pairs.size();
    S.ycross.resize(NP);
    for (int q = 0; q < NP; ++q) S.ycross[q] = S.Y.col(pairs[q].first).dot(S.Y.col(pairs[q].second));
    if (!binary || ctx.braw.empty()) return;
    // exact conditional Bernoulli saddlepoint inputs (see bern_spa_solve). The
    // gap delta cancels identically in the p-value, so the ORIGINAL phenotype's
    // delta is used for every variant; kappa_22 and cylo are variant-specific.
    S.cb.k22.assign(NP, 0.0); S.cb.cdelta.assign(NP, 0.0);
    S.cb.cond_t.assign(NP, 0); S.cb.rand_t.assign(NP, 0);
    S.cb.cpi.assign(NP, std::vector<float>(n, 0.f));
    S.cb.cylo.assign(NP, std::vector<float>(n, 0.f));
    for (int q = 0; q < NP; ++q) {
        const int a = pairs[q].first, b = pairs[q].second;
        double ma = 0, mb = 0;
        for (int i = 0; i < n; ++i) { ma += S.Y(i, a); mb += S.Y(i, b); }
        ma /= n; mb /= n;
        double m22 = 0, m20 = 0, m02 = 0, m11 = 0, p11 = 0;
        for (int i = 0; i < n; ++i) {
            const double x = S.Y(i, a) - ma, y = S.Y(i, b) - mb;
            m22 += x*x*y*y; m20 += x*x; m02 += y*y; m11 += x*y;
            if (ctx.braw[a][i] && ctx.braw[b][i]) p11 += 1.0;
        }
        m22 /= n; m20 /= n; m02 /= n; m11 /= n; p11 /= n;
        S.cb.k22[q] = m22 - m20 * m02 - 2.0 * m11 * m11;
        const int tc = (ctx.prev[a] <= ctx.prev[b]) ? a : b, tr = (tc == a) ? b : a;
        S.cb.cond_t[q] = tc; S.cb.rand_t[q] = tr; S.cb.cdelta[q] = ctx.bdelta[tr];
        const double pc = ctx.prev[tc], prr = ctx.prev[tr];
        const double pi1 = (pc > 0.0) ? p11 / pc : prr;
        const double pi0 = (pc < 1.0) ? (prr - p11) / (1.0 - pc) : prr;
        for (int i = 0; i < n; ++i) {
            double pv = ctx.braw[tc][i] ? pi1 : pi0;
            if (pv < 1e-8) pv = 1e-8; if (pv > 1.0 - 1e-8) pv = 1.0 - 1e-8;
            S.cb.cpi[q][i]  = (float) pv;
            S.cb.cylo[q][i] = (float)(S.Y(i, tr) - (double) ctx.braw[tr][i] * ctx.bdelta[tr]);
        }
    }
    S.cb.active = true;
}


// ===========================================================================
// One assembled unit: [ cat_0 | ... | cat_{A-1} | flank ]
// ===========================================================================
// Used for BOTH passes. `a..b` is the unit's own SNP range; the flanks are two
// explicit index ranges (empty for pass 1; the neighbouring windows in pass 2).
struct ChunkDataA {
    std::string chr; long start, end;
    std::vector<int> cat_m;            // columns per ACTIVE category
    std::vector<int> cat_id;           // global category index of each active category
    std::vector<std::string> cat_name;
    std::vector<double> cat_w;         // sum_j [2f(1-f)]^(1+alpha) per active category:
                                       // the alpha-model's expected variance share, used
                                       // to distribute the window adjustment over chunks
    int m_flank, K;
    GenoMat V;
    std::vector<double> cj;
    ChunkDataA() : m_flank(0), K(0) {}
};
struct ChunkResultA {
    std::string chr; long start, end;
    int m_flank;
    std::vector<std::string> cat_name;
    std::vector<int> cat_m, cat_id;
    std::vector<double> cat_w;
    std::vector< std::vector<double> > vg, se_vg, h2, p_spa;   // [cat][trait or pair]
    std::vector< std::vector<int> > spa_used;
    std::vector<double> vg_flank, vg_env;                       // per trait/pair
    std::vector< std::vector<double> > vg_t1, vg_t2;            // coher only
};

// Read a SNP index range, project covariates, apply alpha; return false on error.
static bool read_block(const ChunkContext& ctx, int lo, int hi, Cell& c) {
    c.X = GenoMat(ctx.n_inds, 0);
    if (hi <= lo) return true;
    if (!read_cell_idx(ctx.bed, ctx.geno_keep, lo, hi, c)) return false;
    if (ctx.covZ.cols() > 0) project_covariates(c.X, ctx.covZ, ctx.covM);
    apply_alpha(c.X, c.maf, ctx.alpha, 0, lo);
    return true;
}

static bool make_chunk(const ChunkContext& ctx, size_t ci, int a, int b,
                       int fL0, int fL1, int fR0, int fR1, ChunkDataA& cd) {
    cd = ChunkDataA();
    cd.chr = ctx.chr_order[ci];
    cd.start = ctx.bim.bp[a]; cd.end = ctx.bim.bp[b - 1];

    Cell tgt;
    if (!read_block(ctx, a, b, tgt)) return false;
    const int m_t = (int) tgt.X.cols();
    if (m_t <= 0) return false;

    const int A = ctx.n_cat;
    std::vector< std::vector<int> > cat_cols(A);
    std::vector<int> bg_cols;
    for (int j = 0; j < m_t; ++j) {
        const int cid = ctx.snp_cat[a + j];
        if (cid >= 0 && cid < A) cat_cols[cid].push_back(j); else bg_cols.push_back(j);
    }
    std::vector<int> act;
    for (int c = 0; c < A; ++c) if (!cat_cols[c].empty()) act.push_back(c);
    if (act.empty()) return false;

    Cell fl, fr;
    if (!read_block(ctx, fL0, fL1, fl)) return false;
    if (!read_block(ctx, fR0, fR1, fr)) return false;
    const int m_flank = (int) bg_cols.size() + (int) fl.X.cols() + (int) fr.X.cols();

    int Kt = 0; for (size_t k = 0; k < act.size(); ++k) Kt += (int) cat_cols[act[k]].size();
    cd.K = Kt + m_flank;
    cd.V = GenoMat(ctx.n_inds, cd.K);
    const double e1 = 1.0 + ctx.alpha;
    int off = 0;
    for (size_t k = 0; k < act.size(); ++k) {
        const int c = act[k]; const std::vector<int>& cc = cat_cols[c];
        double w = 0.0;
        for (size_t q = 0; q < cc.size(); ++q) {
            cd.V.col(off + (int) q) = tgt.X.col(cc[q]);
            const double f = tgt.maf[cc[q]];
            w += (f > 0.0 && f < 1.0) ? std::pow(2.0 * f * (1.0 - f), e1) : 0.0;
        }
        off += (int) cc.size();
        cd.cat_m.push_back((int) cc.size()); cd.cat_id.push_back(c);
        cd.cat_name.push_back(ctx.cat_names[c]); cd.cat_w.push_back(w);
    }
    for (size_t q = 0; q < bg_cols.size(); ++q) cd.V.col(off++) = tgt.X.col(bg_cols[q]);
    if (fl.X.cols() > 0) { cd.V.middleCols(off, fl.X.cols()) = fl.X; off += (int) fl.X.cols(); }
    if (fr.X.cols() > 0) { cd.V.middleCols(off, fr.X.cols()) = fr.X; off += (int) fr.X.cols(); }
    cd.m_flank = m_flank;
    cd.cj.resize(cd.K);
    for (int j = 0; j < cd.K; ++j) cd.cj[j] = cd.V.col(j).cast<double>().squaredNorm();
    return true;
}


static void test_chunk_annot(const ChunkDataA& cd, const Eigen::MatrixXd& Y, const GenoMat& Yf,
                             const std::vector<double>& Vp, const std::vector<double>& yty,
                             const std::vector<double>& kur, bool binary,
                             bool spa, double spa_thresh, int cov_df, bool coher,
                             const std::vector< std::pair<int,int> >& pairs,
                             const std::vector<double>& ycross,
                             const CoherBin& cb,
                             ChunkResultA& cr) {
    const int n = (int) Y.rows(), P = (int) Y.cols(), K = cd.K;
    const int A = (int) cd.cat_m.size();
    const bool has_f = (cd.m_flank > 0);
    // Components: the A tested categories, then the flank ONLY if present. With no
    // flank (pass 1) the model is { cat_0, ..., cat_{A-1}, sigma_e I }.
    const int C = A + (has_f ? 1 : 0), env = C;

    cr.chr = cd.chr; cr.start = cd.start; cr.end = cd.end;
    cr.m_flank = cd.m_flank;
    cr.cat_name = cd.cat_name; cr.cat_m = cd.cat_m; cr.cat_id = cd.cat_id; cr.cat_w = cd.cat_w;
    const int NO = coher ? (int) pairs.size() : P;
    cr.vg_flank.assign(NO, NA_REAL);
    cr.vg_env.assign(NO, NA_REAL);
    cr.vg.assign(A, std::vector<double>(NO, NA_REAL));
    cr.se_vg.assign(A, std::vector<double>(NO, NA_REAL));
    cr.h2.assign(A, std::vector<double>(NO, NA_REAL));
    cr.p_spa.assign(A, std::vector<double>(NO, NA_REAL));
    cr.spa_used.assign(A, std::vector<int>(NO, 0));
    cr.vg_t1.assign(A, std::vector<double>(NO, NA_REAL));
    cr.vg_t2.assign(A, std::vector<double>(NO, NA_REAL));
    if (A <= 0) return;

    // component column offsets in V: cats, then flank (if any), then common (if any)
    std::vector<int> off(C + 1, 0);
    for (int c = 0; c < A; ++c) off[c + 1] = off[c] + cd.cat_m[c];
    int nx = A;
    if (has_f) { off[nx + 1] = off[nx] + cd.m_flank;  ++nx; }

    GenoMat Gf = GenoMat::Zero(K, K);
    Gf.selfadjointView<Eigen::Upper>().rankUpdate(cd.V.transpose());
    Eigen::MatrixXd G = GenoMat(Gf.selfadjointView<Eigen::Upper>()).cast<double>();
    Eigen::MatrixXd G2 = G.array().square();

    std::vector<double> tr(C, 0.0), g(C, 0.0);
    for (int c = 0; c < C; ++c) {
        for (int j = off[c]; j < off[c + 1]; ++j) tr[c] += cd.cj[j];
        if (!(tr[c] > 0.0)) return;
        g[c] = (double) n / tr[c];
    }

    Eigen::MatrixXd T = Eigen::MatrixXd::Zero(C + 1, C + 1);
    for (int a = 0; a < C; ++a)
        for (int b = a; b < C; ++b) {
            double v = g[a] * g[b] * G2.block(off[a], off[b], off[a+1]-off[a], off[b+1]-off[b]).sum();
            T(a, b) = v; T(b, a) = v;
        }
    // See test_chunk: tr(P) = n - cov_df, not n, once covariates are regressed out.
    for (int a = 0; a < C; ++a) { T(a, env) = (double) n; T(env, a) = (double) n; }
    T(env, env) = (double) (n - cov_df);
    Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> Tcod(T);
    Eigen::MatrixXd Tinv = Tcod.pseudoInverse();

    Eigen::MatrixXd L; bool have_L = false;
    if (spa) {
        double gscale = G.trace() / K;
        Eigen::MatrixXd Gj = G; double applied = 0.0;
        for (double r = 1e-4; r <= 1e-1; r *= 10.0) {
            Gj.diagonal().array() += (r - applied) * gscale; applied = r;
            Eigen::LLT<Eigen::MatrixXd> llt(Gj);
            if (llt.info() == Eigen::Success) { L = llt.matrixL(); have_L = true; break; }
        }
    }

    // ---- block outer products, once per chunk (shared across traits AND
    // categories -- the saving is A x P fold here, not just P fold) ----------
    std::vector<Eigen::MatrixXd> W;
    bool use_wcache = spa && have_L && ((double) C * K * K * 8.0 <= 6e8);
    if (use_wcache) {
        W.resize(C);
        for (int a = 0; a < C; ++a) {
            const int ma = off[a + 1] - off[a];
            W[a].noalias() = L.middleRows(off[a], ma).transpose() * L.middleRows(off[a], ma);
        }
    }
    std::vector<double> d0v(C), dmv(C);
    Eigen::MatrixXd Smat, Amat;
    auto build_SA = [&]() {
        if (use_wcache) {
            Smat.setZero(K, K); Amat.setZero(K, K);
            for (int a = 0; a < C; ++a) {
                if (d0v[a] != 0.0) Smat.noalias() += d0v[a] * W[a];
                if (dmv[a] != 0.0) Amat.noalias() += dmv[a] * W[a];
            }
        } else {
            Eigen::VectorXd D0(K), DM(K);
            for (int a = 0; a < C; ++a)
                for (int j = off[a]; j < off[a + 1]; ++j) { D0[j] = d0v[a]; DM[j] = dmv[a]; }
            Smat.noalias() = L.transpose() * (D0.asDiagonal() * L);
            Amat.noalias() = L.transpose() * (DM.asDiagonal() * L);
        }
    };

    // ---- binary: sum_i B_ii^2 per TESTED CATEGORY, once per chunk ----------
    // M_0 differs between categories (row c of Tinv), so one diagonal per
    // category -- but still independent of the trait, so A passes not A*P.
    std::vector<double> Sb2(A, 0.0);
    // Sb2 must be EAGER: it corrects var0, which feeds the Wald screen itself.
    // It is only O(nK), so that is cheap.
    if (binary && spa && have_L) {
        for (int c = 0; c < A; ++c) {
            Eigen::VectorXd bdiag = Eigen::VectorXd::Constant(n, Tinv(c, env));
            for (int a = 0; a < C; ++a) {
                const double dm = Tinv(c, a) * g[a];
                if (dm == 0.0) continue;
                for (int j = off[a]; j < off[a + 1]; ++j)
                    bdiag.array() += dm * cd.V.col(j).cast<double>().array().square();
            }
            Sb2[c] = bdiag.squaredNorm();
        }
    }
    // The 4-cumulant setup is LAZY and PER CATEGORY: an eigensolve with vectors
    // plus an O(nK^2) product, so doing it eagerly for every category of every
    // chunk multiplied the run time by roughly (1 + A). Now a category pays it
    // only when one of its tests actually reaches the saddlepoint.
    // Six cumulants of Q are matched, which needs cumulants of u_k^2 up to order
    // 6, hence moments of u_k up to order 12, hence power sums to order 12.
    const int NPOW = 12;
    const int NCUM = 6;
    std::vector<char> cum4_tried(A, 0), have_cum4(A, 0);
    std::vector<Eigen::VectorXd> lamB(A);
    std::vector<Eigen::MatrixXd> Spow(A);
    auto ensure_cum4 = [&](int c) -> bool {
        if (cum4_tried[c]) return have_cum4[c] != 0;
        cum4_tried[c] = 1;
        if (!(binary && spa && have_L) || kur.empty()) return false;
        {
            Eigen::MatrixXd Qm = Eigen::MatrixXd::Zero(K, K);
            for (int a = 0; a < C; ++a) {
                const double dm = Tinv(c, a) * g[a];
                if (dm != 0.0) Qm.noalias() += dm * W[a];
            }
            Qm = (0.5 * (Qm + Qm.transpose())).eval();
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esQ(Qm);
            if (esQ.info() != Eigen::Success) return false;
            lamB[c] = esQ.eigenvalues().array() + Tinv(c, env);
            Eigen::MatrixXd Wmat =
                L.transpose().triangularView<Eigen::Upper>().solve(esQ.eigenvectors());
            Spow[c] = Eigen::MatrixXd::Zero(K, NPOW + 1);
            Eigen::VectorXd nrm2 = Eigen::VectorXd::Zero(K);
            const int BS = 4096;
            for (int r0 = 0; r0 < n; r0 += BS) {
                const int rr = std::min(BS, n - r0);
                Eigen::MatrixXd Ablk = cd.V.middleRows(r0, rr).cast<double>() * Wmat;
                nrm2 += Ablk.colwise().squaredNorm();
                Eigen::MatrixXd Pw = Ablk;
                for (int r = 2; r <= NPOW; ++r) {
                    Pw = Pw.cwiseProduct(Ablk);
                    Spow[c].col(r) += Pw.colwise().sum().transpose();
                }
            }
            bool okc = true;
            for (int k = 0; k < K; ++k) {
                const double nn2 = nrm2[k];
                if (!(nn2 > 0.0) || !std::isfinite(nn2)) { okc = false; break; }
                const double nn = std::sqrt(nn2); double pw = nn2;
                for (int r = 2; r <= NPOW; ++r) { Spow[c](k, r) /= pw; pw *= nn; }
            }
            have_cum4[c] = okc ? 1 : 0;
        }
        return have_cum4[c] != 0;
    };

    Eigen::MatrixXd U = (cd.V.transpose() * Yf).cast<double>();      // K x P, one GEMM

    // ---- CO-HERITABILITY BRANCH (annotation path) --------------------------
    // Same construction as the plain path, but every one of the A tested
    // categories has its own M_c, hence its own A_c and its own null spectrum.
    if (coher) {
        std::vector<Eigen::VectorXd> sg(P);
        for (int t = 0; t < P; ++t) {
            Eigen::VectorXd qu(C + 1);
            for (int cc = 0; cc < C; ++cc) {
                double qa = 0.0;
                for (int j = off[cc]; j < off[cc + 1]; ++j) qa += U(j, t) * U(j, t);
                qu[cc] = g[cc] * qa;
            }
            qu[env] = yty[t];
            sg[t] = Tcod.solve(qu);
        }
        std::vector<Eigen::MatrixXd> Acat(A);
        if (spa && have_L) {
            for (int c = 0; c < A; ++c) {
                std::vector<double> dm(C);
                for (int a2 = 0; a2 < C; ++a2) dm[a2] = Tinv(c, a2) * g[a2];
                coher_block(L, W, use_wcache, off, C, K, dm, Tinv(c, env), Acat[c]);
            }
        }
        const bool spec_ok = spa && have_L &&
            ((double) 4 * K * K * 8.0 * 3.0 <= 2.0e9);

        for (int pi = 0; pi < (int) pairs.size(); ++pi) {
            const int t1 = pairs[pi].first, t2 = pairs[pi].second;
            Eigen::VectorXd q(C + 1);
            for (int cc = 0; cc < C; ++cc) {
                double qa = 0.0;
                for (int j = off[cc]; j < off[cc + 1]; ++j) qa += U(j, t1) * U(j, t2);
                q[cc] = g[cc] * qa;
            }
            q[env] = ycross[pi];
            Eigen::VectorXd s12 = Tcod.solve(q);

            const double vv = Vp[t1] * Vp[t2];
            for (int c = 0; c < A; ++c) {
                cr.vg[c][pi]    = s12[c];
                cr.vg_t1[c][pi] = sg[t1][c];
                cr.vg_t2[c][pi] = sg[t2][c];
                cr.h2[c][pi] = (vv > 0.0) ? s12[c] / std::sqrt(vv) : NA_REAL;
            }
            cr.vg_flank[pi] = has_f ? s12[A] : NA_REAL;
            cr.vg_env[pi]   = s12[env];
            if (!spa || !have_L) continue;

            for (int c = 0; c < A; ++c) {
                CoherFit F;
                double Sb2c = 0.0;
                if (binary && cb.active) {
                    Eigen::VectorXd md; coher_mdiag(cd.V, off, C, g, Tinv, c, env, n, md);
                    Sb2c = md.squaredNorm();
                }
                coher_fit(c, Tinv, C, env, g, off, K, n, L, W, use_wcache, Acat[c],
                          sg[t1], sg[t2], s12, Vp[t1], Vp[t2],
                          (binary && cb.active) ? cb.k22[pi] : 0.0, Sb2c, F);
                if (!F.ok) continue;
                cr.se_vg[c][pi] = std::sqrt(F.var0);
                const double p_wald =
                    std::erfc(std::abs(s12[c] / cr.se_vg[c][pi]) / std::sqrt(2.0));
                // ---- BINARY: exact conditional saddlepoint ------------------
                if (binary && cb.active) {
                    const int tc = cb.cond_t[pi], trd = cb.rand_t[pi];
                    Eigen::VectorXd wv;
                    coher_wvec(cd.V, off, C, K, g, Tinv, c, env, U.col(tc),
                                   Y.col(tc), wv);
                    std::vector<double> cvec(n);
                    double A = 0.0;
                    const double dl = cb.cdelta[pi];
                    for (int i = 0; i < n; ++i) {
                        A += wv[i] * cb.cylo[pi][i];
                        cvec[i] = dl * wv[i];
                    }
                    QuadSpaResult br = bern_spa_solve(s12[c], cvec, cb.cpi[pi], A);
                    if (br.converged) { cr.p_spa[c][pi] = br.p; cr.spa_used[c][pi] = 4; }
                    else                  { cr.p_spa[c][pi] = p_wald; cr.spa_used[c][pi] = 0; }
                    continue;
                }
                if (p_wald >= spa_thresh || !spec_ok) {
                    cr.p_spa[c][pi] = p_wald; cr.spa_used[c][pi] = 0; continue;
                }
                std::vector<double> eig; double erep, l2, nrep;
                if (!coher_spectrum(F, K, n, eig, erep, l2, nrep)) {
                    cr.p_spa[c][pi] = p_wald; cr.spa_used[c][pi] = 0; continue;
                }
                QuadSpaResult qr = quad_spa_solve(s12[c], eig, erep, nrep,
                                                  100, 1e-8, l2, nrep);
                if (qr.converged) { cr.p_spa[c][pi] = qr.p; cr.spa_used[c][pi] = 1; }
                else              { cr.p_spa[c][pi] = p_wald; cr.spa_used[c][pi] = 0; }
            }
        }
        return;
    }


    for (int t = 0; t < P; ++t) {
        const Eigen::VectorXd u = U.col(t);
        Eigen::VectorXd q(C + 1);
        for (int c = 0; c < C; ++c) {
            double qc = 0.0; for (int j = off[c]; j < off[c + 1]; ++j) qc += u[j] * u[j];
            q[c] = g[c] * qc;
        }
        q[env] = yty[t];
        Eigen::VectorXd sigma = Tcod.solve(q);

        for (int c = 0; c < A; ++c) {
            cr.vg[c][t] = sigma[c];
            cr.h2[c][t] = (Vp[t] > 0) ? sigma[c] / Vp[t] : NA_REAL;
        }
        cr.vg_flank[t] = has_f ? sigma[A] : NA_REAL;
        cr.vg_env[t]   = sigma[env];
        if (!spa || !have_L) continue;

        double sigma_env0 = std::max(sigma[env], 1e-8 * (Vp[t] > 0 ? Vp[t] : 1.0));
        Eigen::VectorXd s0base = sigma.head(C);
        for (int a = 0; a < C; ++a) if (s0base[a] < 0.0) s0base[a] = 0.0;

        // ---- SPA-test each category component c ----
        for (int c = 0; c < A; ++c) {
            const double c_env = Tinv(c, env);
            for (int a = 0; a < C; ++a) {
                d0v[a] = (a == c) ? 0.0 : s0base[a] * g[a];   // null the tested category
                dmv[a] = Tinv(c, a) * g[a];
            }
            build_SA();
            Smat.diagonal().array() += sigma_env0;
            Amat.diagonal().array() += c_env;

            // Cholesky-free Wald screen: 2 tr((A S)^2) = 2 ||Cf^T A Cf||_F^2
            double eig_rep = c_env * sigma_env0;
            long n_rep = (long) n - K;
            Eigen::MatrixXd N; N.noalias() = Amat * Smat;
            double var0 = 2.0 * (N.cwiseProduct(N.transpose())).sum();
            if (n_rep > 0) var0 += (double) n_rep * 2.0 * eig_rep * eig_rep;
            if (!(var0 > 0.0)) continue;

            double cscale = 1.0;                       // binary variance correction
            if (binary && Sb2[c] > 0.0) {
                const double d2 = kur[(size_t) t * (MAXCUM + 1) + 4] * Sb2[c];
                if (var0 + d2 > 0.0) { cscale = std::sqrt((var0 + d2) / var0); var0 += d2; }
            }

            cr.se_vg[c][t] = std::sqrt(var0);
            double p_wald = std::erfc(std::abs(sigma[c] / cr.se_vg[c][t]) / std::sqrt(2.0));
            if (p_wald >= spa_thresh) { cr.p_spa[c][t] = p_wald; cr.spa_used[c][t] = 0; continue; }

            // 4-cumulant path first: when it succeeds, none of the Cholesky of
            // Smat, the Cf products, or the eigensolve below is needed.
            bool used_cum4 = false;
            if (binary && ensure_cum4(c) && lamB[c].size() == K) {
                const double* kz = &kur[(size_t) t * (MAXCUM + 1)];
                double KQ[NCUM + 1] = {0,0,0,0,0,0,0};
                bool ok4 = true;
                for (int k = 0; k < K && ok4; ++k) {
                    double ku[MAXCUM + 1] = {0};
                    for (int r = 2; r <= NPOW; ++r) ku[r] = kz[r] * Spow[c](k, r);
                    double mu[MAXCUM + 1]; cums_to_moms(ku, mu, NPOW);
                    double ms[NCUM + 1]; ms[0] = 1.0;
                    for (int s2 = 1; s2 <= NCUM; ++s2) ms[s2] = mu[2 * s2];
                    double csq[NCUM + 1]; moms_to_cums(ms, csq, NCUM);
                    double lp = 1.0;
                    for (int mm = 1; mm <= NCUM; ++mm) {
                        lp *= lamB[c][k];
                        if (!std::isfinite(csq[mm])) { ok4 = false; break; }
                        KQ[mm] += lp * csq[mm];
                    }
                }
                if (ok4) {
                    const double cev = Tinv(c, env), Vpt = (Vp[t] > 0 ? Vp[t] : 1.0);
                    const double fac[NCUM + 1] = {0, 1, 2, 8, 48, 384, 3840};
                    double lp = 1.0, vp = 1.0;
                    for (int mm = 1; mm <= NCUM; ++mm) {
                        lp *= cev; vp *= Vpt;
                        KQ[mm] += fac[mm] * (double)(n - K) * lp * vp;
                    }
                    // Shape only -- pin kappa2 to the exact var0. See test_chunk.
                    if (KQ[2] > 0.0 && var0 > 0.0 && std::isfinite(KQ[2])) {
                        const double s = std::sqrt(var0 / KQ[2]);
                        double sp = 1.0;
                        for (int mm = 1; mm <= NCUM; ++mm) { sp *= s; KQ[mm] *= sp; }
                    }
                    if (std::isfinite(KQ[2]) && KQ[2] > 0) {
                        std::vector<double> none;
                        double L3[3], N3[3];
                        if (fit_three_component(KQ, L3, N3)) {
                            QuadSpaResult q6 = quad_spa_solve(sigma[c], none, L3[0], N3[0],
                                                              100, 1e-8, L3[1], N3[1],
                                                              L3[2], N3[2]);
                            if (q6.converged) {
                                cr.p_spa[c][t] = q6.p; cr.spa_used[c][t] = 3; used_cum4 = true;
                            }
                        }
                        double l1, v1, l2, v2;
                        if (!used_cum4 && fit_two_component(KQ, &l1, &v1, &l2, &v2)) {
                            QuadSpaResult q4 = quad_spa_solve(sigma[c], none, l1, v1,
                                                              100, 1e-8, l2, v2);
                            if (q4.converged) {
                                cr.p_spa[c][t] = q4.p; cr.spa_used[c][t] = 2; used_cum4 = true;
                            }
                        }
                    }
                }
            }
            if (used_cum4) continue;

            // ---- fallback: build Asym and solve the eigenproblem -------------
            Eigen::LLT<Eigen::MatrixXd> lltS(Smat);
            if (lltS.info() != Eigen::Success) { cr.p_spa[c][t] = p_wald; continue; }
            Eigen::MatrixXd Cf = lltS.matrixL();
            Eigen::MatrixXd Asym = Cf.transpose() * Amat * Cf;
            Asym = (0.5 * (Asym + Asym.transpose())).eval();
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ses(Asym, Eigen::EigenvaluesOnly);
            if (ses.info() != Eigen::Success) { cr.p_spa[c][t] = p_wald; continue; }
            const Eigen::VectorXd& ev = ses.eigenvalues();
            std::vector<double> eig(ev.data(), ev.data() + K);
            double eig_rep_s = eig_rep;
            if (cscale != 1.0) {
                for (size_t e = 0; e < eig.size(); ++e) eig[e] *= cscale;
                eig_rep_s *= cscale;
            }
            QuadSpaResult qr = quad_spa_solve(sigma[c], eig, eig_rep_s, n_rep);
            if (qr.converged) { cr.p_spa[c][t] = qr.p; cr.spa_used[c][t] = 1; }
            else              { cr.p_spa[c][t] = p_wald; cr.spa_used[c][t] = 0; }
        }
    }
}


// ===========================================================================
// Driver
// ===========================================================================
struct ChunkParams {
    int chunk_size; long window_bp; bool do_windows;
    bool spa; double spa_thresh; bool binary;
    int  cov_df;
    std::string out_file; int batch_size, n_threads;
    std::string chr;
    bool coher;
    std::vector< std::pair<int,int> > pairs;
    std::vector<std::string> pair_names;
};

// One unit of work: its own range [a,b) plus two flank ranges. `win` is the
// index of the window a pass-1 chunk belongs to (or of the window itself).
struct Job { size_t ci; int a, b, fL0, fL1, fR0, fR1; int win; };

struct Worker : public RcppParallel::Worker {
    const ChunkContext& ctx; const std::vector<Job>& jobs; size_t job0;
    const PhenoStats& S; const ChunkParams& pr; bool spa;
    std::vector<ChunkResultA>& out; std::vector<char>& ok;
    Worker(const ChunkContext& ctx, const std::vector<Job>& jobs, size_t job0,
           const PhenoStats& S, const ChunkParams& pr, bool spa,
           std::vector<ChunkResultA>& out, std::vector<char>& ok)
        : ctx(ctx), jobs(jobs), job0(job0), S(S), pr(pr), spa(spa), out(out), ok(ok) {}
    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t w = begin; w < end; ++w) {
            const Job& j = jobs[job0 + w];
            ChunkDataA cd;
            if (!make_chunk(ctx, j.ci, j.a, j.b, j.fL0, j.fL1, j.fR0, j.fR1, cd)) { ok[w] = 0; continue; }
            test_chunk_annot(cd, S.Y, S.Yf, S.Vp, S.yty, S.kur, pr.binary, spa, pr.spa_thresh,
                             pr.cov_df, pr.coher, pr.pairs, S.ycross, S.cb, out[w]);
            ok[w] = 1;
        }
    }
};

static void run_jobs(const ChunkContext& ctx, const ChunkParams& pr, const PhenoStats& S,
                     bool spa, const std::vector<Job>& jobs, size_t j0, size_t j1,
                     std::vector<ChunkResultA>& res, std::vector<char>& ok,
                     const char* tag, long& done_total, long total) {
    typedef std::chrono::steady_clock clk;
    static clk::time_point t0 = clk::now();
    size_t idx = j0;
    while (idx < j1) {
        const size_t nb = std::min((size_t) std::max(1, pr.batch_size), j1 - idx);
        std::vector<ChunkResultA> tmp(nb); std::vector<char> tok(nb, 0);
        Worker wk(ctx, jobs, idx, S, pr, spa, tmp, tok);
        RcppParallel::parallelFor(0, nb, wk);
        for (size_t b = 0; b < nb; ++b) { res[idx + b] = tmp[b]; ok[idx + b] = tok[b]; }
        idx += nb; done_total += (long) nb;
        double el = std::chrono::duration<double>(clk::now() - t0).count();
        double eta = (done_total > 0) ? el * ((double) total - done_total) / done_total : 0.0;
        Rcout << "\r[" << tag << "] " << done_total << "/" << total << "  "
              << (int)(100.0 * done_total / std::max(1L, total)) << "%  elapsed " << (int) el
              << "s  eta " << (int) eta << "s     " << std::flush;
        Rcpp::checkUserInterrupt();
    }
}

static void report_plan(const ChunkContext& ctx, const ChunkParams& pr, int Kmax, const char* what) {
    int th = pr.n_threads;
    if (th <= 0) { const char* e = std::getenv("RCPP_PARALLEL_NUM_THREADS"); th = e ? std::atoi(e) : 1; if (th <= 0) th = 1; }
    const double V_mb = (double) ctx.n_inds * Kmax * 4.0 / 1048576.0;
    const double K_mb = (double) Kmax * Kmax * 8.0 / 1048576.0;
    const double per  = V_mb + K_mb * (ctx.n_cat + 6);
    Rcout << "  " << what << ": K up to ~" << Kmax << ", ~" << (int) per << " MB per thread, ~"
          << (per * th / 1024.0) << " GB at " << th << " threads.\n";
}

static Rcpp::List run_all(ChunkContext& ctx, const ChunkParams& pr) {
    const int P = ctx.n_pheno, n = ctx.n_inds;
    Eigen::setNbThreads(1);
    if (pr.n_threads > 0) {
        std::string nt = std::to_string(pr.n_threads);
#ifdef _WIN32
        _putenv_s("RCPP_PARALLEL_NUM_THREADS", nt.c_str());
#else
        setenv("RCPP_PARALLEL_NUM_THREADS", nt.c_str(), 1);
#endif
    }
    // h2 denominators: the ORIGINAL phenotypic variance, so pass-1 and pass-2
    // heritabilities (and the adjusted ones) are all fractions of the same thing.
    std::vector<double> Vp0(P);
    for (int t = 0; t < P; ++t) { double m = ctx.Y.col(t).mean();
        Vp0[t] = (ctx.Y.col(t).array() - m).square().sum() / (n - 1); }
    const int NOUT = pr.coher ? (int) pr.pairs.size() : P;
    auto label = [&](int t) -> std::string {
        return pr.coher ? pr.pair_names[t] : as<std::string>(ctx.trait_names[t]); };
    auto h2den = [&](int t) -> double {
        if (!pr.coher) return Vp0[t];
        return std::sqrt(Vp0[pr.pairs[t].first] * Vp0[pr.pairs[t].second]); };

    // ---- enumerate chunks and windows ------------------------------------
    std::vector<Job> cj, wj;                    // chunk jobs, window jobs
    std::vector<size_t> chr_c0, chr_c1;          // chunk job range per chromosome
    std::vector<int> win_ci;                     // chromosome of each window
    for (size_t ci = 0; ci < ctx.chr_order.size(); ++ci) {
        if (!pr.chr.empty() && strip_chr(ctx.chr_order[ci]) != strip_chr(pr.chr)) continue;
        const int lo = ctx.chr_lo[ci], hi = ctx.chr_hi[ci];
        chr_c0.push_back(cj.size());
        // windows: consecutive whole chunks until the bp span reaches window_bp
        std::vector< std::pair<int,int> > wr;    // [a, b) of each window
        int wa = lo;
        for (int a = lo; a <= hi; ) {
            const int b = std::min(hi + 1, a + pr.chunk_size);
            Job j; j.ci = ci; j.a = a; j.b = b; j.fL0 = j.fL1 = j.fR0 = j.fR1 = 0;
            j.win = (int) wj.size() + (int) wr.size();
            cj.push_back(j);
            const bool last = (b > hi);
            if (last || ctx.bim.bp[b - 1] - ctx.bim.bp[wa] + 1 >= pr.window_bp) {
                wr.push_back(std::make_pair(wa, b)); wa = b;
            }
            a = b;
        }
        // A trailing remainder shorter than half a window is absorbed into its
        // predecessor (as the old build_cells did), so no window-level estimate
        // rests on a sliver. The chunks it held are re-pointed accordingly.
        if (wr.size() >= 2) {
            const std::pair<int,int>& lastw = wr.back();
            if (ctx.bim.bp[lastw.second - 1] - ctx.bim.bp[lastw.first] + 1 < pr.window_bp / 2) {
                const int wlast = (int) wj.size() + (int) wr.size() - 1;
                wr[wr.size() - 2].second = lastw.second; wr.pop_back();
                for (size_t k = chr_c0.back(); k < cj.size(); ++k) if (cj[k].win == wlast) cj[k].win = wlast - 1;
            }
        }
        for (size_t w = 0; w < wr.size(); ++w) {
            Job j; j.ci = ci; j.a = wr[w].first; j.b = wr[w].second;
            j.fL0 = (w > 0) ? wr[w - 1].first : j.a;               j.fL1 = j.a;
            j.fR0 = j.b; j.fR1 = (w + 1 < wr.size()) ? wr[w + 1].second : j.b;
            j.win = (int) wj.size(); wj.push_back(j); win_ci.push_back((int) ci);
        }
        chr_c1.push_back(cj.size());
    }
    Rcout << "Pass 1: " << cj.size() << " chunks of " << pr.chunk_size << " SNPs, no flanks, "
          << ctx.n_cat << " categories, phenotype = y - LOCO PRS of the tested chromosome.\n";
    report_plan(ctx, pr, pr.chunk_size, "pass 1");
    if (pr.do_windows) {
        int Kmax = 0; for (size_t w = 0; w < wj.size(); ++w)
            Kmax = std::max(Kmax, (wj[w].b - wj[w].a) + (wj[w].fL1 - wj[w].fL0) + (wj[w].fR1 - wj[w].fR0));
        Rcout << "Pass 2: " << wj.size() << " windows of ~" << pr.window_bp
              << " bp with neighbouring windows as flanks, phenotype = y - full PRS per prs_mask.\n";
        report_plan(ctx, pr, Kmax, "pass 2");
    }

    // ---- PASS 1: association, one phenotype variant per chromosome ---------
    std::vector<ChunkResultA> r1(cj.size()); std::vector<char> ok1(cj.size(), 0);
    long done = 0;
    for (size_t k = 0; k < chr_c0.size(); ++k) {
        if (chr_c1[k] == chr_c0[k]) continue;
        const size_t ci = cj[chr_c0[k]].ci;
        Eigen::MatrixXd off = Eigen::MatrixXd::Zero(n, P);
        // ctx.prs is stored trait-major: prs[t * G + g]
        const size_t G = ctx.prs_mask.empty() ? 0 : ctx.prs_mask[0].size();
        for (int t = 0; t < P; ++t)
            for (size_t g = 0; g < G; ++g) {
                const PrsSource& S = ctx.prs[(size_t) t * G + g];
                int col = -1;
                for (size_t c = 0; c < S.chr_index.size(); ++c) if (S.chr_index[c] == (int) ci) { col = (int) c; break; }
                if (col >= 0) off.col(t) += S.loco.col(col);
                else          off.col(t) += S.full;     // chromosome absent from the fit
            }
        PhenoStats S; build_pheno_stats(ctx, off, pr.binary, pr.coher, pr.pairs, S);
        run_jobs(ctx, pr, S, pr.spa, cj, chr_c0[k], chr_c1[k], r1, ok1, "pass1", done, (long) cj.size());
    }
    Rcout << "\n";

    // ---- PASS 2: windows, one phenotype variant per distinct prs_mask row ---
    const size_t G = ctx.prs_mask.empty() ? 0 : ctx.prs_mask[0].size();
    std::vector<int> cat_variant(ctx.n_cat, 0);
    std::vector< std::vector<unsigned char> > variants;           // distinct mask rows
    for (int c = 0; c < ctx.n_cat; ++c) {
        const std::vector<unsigned char>& row = ctx.prs_mask[c];
        int v = -1;
        for (size_t q = 0; q < variants.size(); ++q) if (variants[q] == row) { v = (int) q; break; }
        if (v < 0) { v = (int) variants.size(); variants.push_back(row); }
        cat_variant[c] = v;
    }
    std::vector< std::vector<ChunkResultA> > r2(variants.size(), std::vector<ChunkResultA>(wj.size()));
    std::vector< std::vector<char> > ok2(variants.size(), std::vector<char>(wj.size(), 0));
    if (pr.do_windows && !wj.empty()) {
        long done2 = 0; const long tot2 = (long) (wj.size() * variants.size());
        for (size_t v = 0; v < variants.size(); ++v) {
            Eigen::MatrixXd off = Eigen::MatrixXd::Zero(n, P);
            for (int t = 0; t < P; ++t)
                for (size_t g = 0; g < G; ++g)
                    if (variants[v][g]) off.col(t) += ctx.prs[(size_t) t * G + g].full;
            PhenoStats S; build_pheno_stats(ctx, off, false, pr.coher, pr.pairs, S);
            run_jobs(ctx, pr, S, false, wj, 0, wj.size(), r2[v], ok2[v], "pass2", done2, tot2);
        }
        Rcout << "\n";
    }

    // ---- ADJUSTMENT: distribute (window - sum of chunks) over the chunks -----
    // adj[chunk][k][t] for active category k of that chunk; NA if no window value
    std::vector< std::vector< std::vector<double> > > adj(cj.size());
    for (size_t k = 0; k < cj.size(); ++k)
        if (ok1[k]) adj[k].assign(r1[k].cat_m.size(), std::vector<double>(NOUT, NA_REAL));
    if (pr.do_windows) {
        // chunks of each window
        std::vector< std::vector<size_t> > members(wj.size());
        for (size_t k = 0; k < cj.size(); ++k) if (ok1[k]) members[cj[k].win].push_back(k);
        for (size_t w = 0; w < wj.size(); ++w) {
            for (int c = 0; c < ctx.n_cat; ++c) {
                const int v = cat_variant[c];
                if (!ok2[v][w]) continue;
                const ChunkResultA& W = r2[v][w];
                int kw = -1; for (size_t q = 0; q < W.cat_id.size(); ++q) if (W.cat_id[q] == c) { kw = (int) q; break; }
                if (kw < 0) continue;
                // sum over member chunks where c is active
                double wsum = 0.0;
                std::vector<double> ssum(NOUT, 0.0);
                std::vector< std::pair<size_t,int> > hits;
                for (size_t m = 0; m < members[w].size(); ++m) {
                    const size_t k = members[w][m]; const ChunkResultA& R = r1[k];
                    for (size_t q = 0; q < R.cat_id.size(); ++q) if (R.cat_id[q] == c) {
                        hits.push_back(std::make_pair(k, (int) q)); wsum += R.cat_w[q];
                        for (int t = 0; t < NOUT; ++t) if (!std::isnan(R.vg[q][t])) ssum[t] += R.vg[q][t];
                    }
                }
                if (hits.empty() || !(wsum > 0.0)) continue;
                for (size_t h = 0; h < hits.size(); ++h) {
                    const size_t k = hits[h].first; const int q = hits[h].second;
                    const double share = r1[k].cat_w[q] / wsum;
                    for (int t = 0; t < NOUT; ++t) {
                        const double Wv = W.vg[kw][t];
                        if (std::isnan(Wv) || std::isnan(r1[k].vg[q][t])) continue;
                        adj[k][q][t] = r1[k].vg[q][t] + (Wv - ssum[t]) * share;
                    }
                }
            }
        }
    }

    // ---- WRITE ---------------------------------------------------------------
    auto wr = [](std::ofstream& f, double v) { if (std::isnan(v)) f << "NA"; else f << v; };
    long n1 = 0, n2 = 0;
    if (!pr.out_file.empty()) {
        std::ofstream f1(pr.out_file.c_str());
        if (!f1.is_open()) stop("Could not open out_file: " + pr.out_file);
        f1 << "chr\tstart\tend\twindow\tcategory\tm_cat\tm_flank\tphenotype\tvg\tse_vg\th2\tvg_flank\tvg_env";
        if (pr.spa)   f1 << "\tp_spa\tspa_used";
        if (pr.coher) f1 << "\tvg_t1\tvg_t2";
        if (pr.do_windows) f1 << "\tvg_adj\th2_adj";
        f1 << "\n";
        for (size_t k = 0; k < cj.size(); ++k) {
            if (!ok1[k]) continue; ++n1;
            const ChunkResultA& R = r1[k];
            for (size_t q = 0; q < R.cat_m.size(); ++q)
                for (int t = 0; t < NOUT; ++t) {
                    f1 << R.chr << '\t' << R.start << '\t' << R.end << '\t' << cj[k].win << '\t'
                       << R.cat_name[q] << '\t' << R.cat_m[q] << '\t' << R.m_flank << '\t' << label(t) << '\t';
                    wr(f1, R.vg[q][t]); f1 << '\t'; wr(f1, R.se_vg[q][t]); f1 << '\t';
                    wr(f1, h2den(t) > 0 ? R.vg[q][t] / h2den(t) : NA_REAL);
                    f1 << '\t'; wr(f1, R.vg_flank[t]); f1 << '\t'; wr(f1, R.vg_env[t]);
                    if (pr.spa) { f1 << '\t'; wr(f1, R.p_spa[q][t]); f1 << '\t' << R.spa_used[q][t]; }
                    if (pr.coher) { f1 << '\t'; wr(f1, R.vg_t1[q][t]); f1 << '\t'; wr(f1, R.vg_t2[q][t]); }
                    if (pr.do_windows) {
                        const double a = adj[k][q][t];
                        f1 << '\t'; wr(f1, a); f1 << '\t'; wr(f1, (!std::isnan(a) && h2den(t) > 0) ? a / h2den(t) : NA_REAL);
                    }
                    f1 << '\n';
                }
        }
        f1.close();
        if (pr.do_windows) {
            std::ofstream f2((pr.out_file + ".windows").c_str());
            if (!f2.is_open()) stop("Could not open " + pr.out_file + ".windows");
            f2 << "chr\tstart\tend\twindow\tn_chunks\tcategory\tm_cat\tm_flank\tprs_variant\tphenotype\tvg\th2\tvg_flank\tvg_env\n";
            std::vector<int> nchunks(wj.size(), 0);
            for (size_t k = 0; k < cj.size(); ++k) if (ok1[k]) ++nchunks[cj[k].win];
            for (size_t w = 0; w < wj.size(); ++w)
                for (int c = 0; c < ctx.n_cat; ++c) {
                    const int v = cat_variant[c]; if (!ok2[v][w]) continue;
                    const ChunkResultA& W = r2[v][w];
                    int kw = -1; for (size_t q = 0; q < W.cat_id.size(); ++q) if (W.cat_id[q] == c) { kw = (int) q; break; }
                    if (kw < 0) continue;
                    ++n2;
                    std::string vs; for (size_t g = 0; g < G; ++g) vs += (variants[v][g] ? '1' : '0');
                    if (vs.empty()) vs = "none";
                    for (int t = 0; t < NOUT; ++t) {
                        f2 << W.chr << '\t' << W.start << '\t' << W.end << '\t' << w << '\t' << nchunks[w] << '\t'
                           << W.cat_name[kw] << '\t' << W.cat_m[kw] << '\t' << W.m_flank << '\t' << vs << '\t'
                           << label(t) << '\t';
                        wr(f2, W.vg[kw][t]); f2 << '\t';
                        wr(f2, h2den(t) > 0 ? W.vg[kw][t] / h2den(t) : NA_REAL);
                        f2 << '\t'; wr(f2, W.vg_flank[t]); f2 << '\t'; wr(f2, W.vg_env[t]); f2 << '\n';
                    }
                }
            f2.close();
        }
    }
    long skipped = 0; for (size_t k = 0; k < cj.size(); ++k) if (!ok1[k]) ++skipped;
    Rcout << "Chunks tested: " << n1 << (skipped ? "  (skipped " + std::to_string(skipped) + ")" : "") << "\n";
    if (pr.do_windows) Rcout << "Window x category estimates: " << n2 << "\n";
    return List::create(_["n_chunks"] = (double) n1, _["n_windows"] = (double) wj.size(),
                        _["trait_names"] = ctx.trait_names, _["categories"] = wrap(ctx.cat_names));
}

}  // end anonymous namespace


// ===========================================================================
// Export
// ===========================================================================
// [[Rcpp::export]]
Rcpp::List stratgwas_run(const std::string& filename,
                         const SEXP pheno_mat,
                         Rcpp::Nullable<Rcpp::IntegerMatrix> annotation = R_NilValue,
                         Rcpp::Nullable<Rcpp::CharacterVector> annot_names = R_NilValue,
                         Rcpp::Nullable<Rcpp::CharacterMatrix> loco_prs = R_NilValue,
                         Rcpp::Nullable<Rcpp::IntegerMatrix> prs_mask = R_NilValue,
                         int chunk_size = 256,
                         double window_bp = 1e6,
                         bool do_windows = true,
                         double alpha = -1.0,
                         Rcpp::Nullable<Rcpp::NumericMatrix> covariates = R_NilValue,
                         double cov_df = NA_REAL,
                         bool SPA = true,
                         double spa_pval_threshold = 0.1,
                         bool binary = false,
                         Rcpp::Nullable<Rcpp::NumericMatrix> binary_raw = R_NilValue,
                         bool coher = false,
                         SEXP chr = R_NilValue,
                         std::string out_file = "",
                         int batch_size = 64,
                         int n_threads = 0) {
    if (chunk_size < 1) stop("chunk_size must be >= 1");
    if (window_bp < chunk_size) stop("window_bp is smaller than one chunk");
    ChunkContext ctx = setup_context(filename, pheno_mat, alpha, covariates, annotation, annot_names);
    const int P = ctx.n_pheno;

    ChunkParams pr;
    pr.chunk_size = chunk_size; pr.window_bp = (long) window_bp; pr.do_windows = do_windows;
    pr.spa = SPA; pr.spa_thresh = spa_pval_threshold; pr.binary = binary;
    pr.cov_df = ISNAN(cov_df) ? ((int) ctx.covZ.cols() + 1) : (int) cov_df;
    if (pr.cov_df < 0) pr.cov_df = 0;
    if (pr.cov_df >= ctx.n_inds) stop("cov_df must be smaller than the sample size");
    Rcout << "Environment moment: T(env,env) = n - " << pr.cov_df << "\n";
    pr.out_file = out_file; pr.batch_size = batch_size; pr.n_threads = n_threads;

    // ---- LOCO PRS sources: a character matrix, P rows (traits) x G columns ----
    // (sources). A single row is recycled over traits with a warning, since a
    // KVIK step-1 fit is per phenotype and sharing one PRS across traits is
    // almost never what you want.
    size_t G = 0;
    if (loco_prs.isNotNull()) {
        Rcpp::CharacterMatrix pm(loco_prs.get());
        G = (size_t) pm.ncol();
        if (pm.nrow() != P && pm.nrow() != 1)
            stop("loco_prs must have one row per phenotype (or a single row to recycle)");
        if (pm.nrow() == 1 && P > 1)
            Rcpp::warning("loco_prs has one row: the same PRS is subtracted from every phenotype");
        ctx.prs.resize((size_t) P * G);
        for (int t = 0; t < P; ++t)
            for (size_t g = 0; g < G; ++g) {
                const int r = (pm.nrow() == 1) ? 0 : t;
                prs_load(as<std::string>(pm(r, (int) g)), ctx, ctx.prs[(size_t) t * G + g]);
            }
    }
    // prs_mask: categories x sources, 1 = subtract that source's FULL PRS when
    // estimating that category in pass 2. Default: subtract everything from
    // everything. Pass 1 always subtracts every source's LOCO PRS.
    ctx.prs_mask.assign(ctx.n_cat, std::vector<unsigned char>(G, 1));
    if (prs_mask.isNotNull()) {
        Rcpp::IntegerMatrix mm(prs_mask.get());
        if (mm.nrow() != ctx.n_cat || (size_t) mm.ncol() != G)
            stop("prs_mask must be n_categories x n_prs_sources");
        for (int c = 0; c < ctx.n_cat; ++c)
            for (size_t g = 0; g < G; ++g) ctx.prs_mask[c][g] = mm(c, (int) g) ? 1 : 0;
    }
    if (G > 0) {
        Rcout << "PRS mask (pass 2, category x source):\n";
        for (int c = 0; c < ctx.n_cat; ++c) {
            Rcout << "  " << ctx.cat_names[c] << ":";
            for (size_t g = 0; g < G; ++g) Rcout << " " << (int) ctx.prs_mask[c][g];
            Rcout << "\n";
        }
    } else if (do_windows) {
        Rcout << "No loco_prs: pass 2 runs on the raw phenotype (no PRS offset).\n";
    }

    // ---- co-heritability: all pairs --------------------------------------------
    pr.coher = coher;
    if (coher) {
        if (P < 2) stop("coher = TRUE needs at least 2 phenotype columns");
        for (int a = 0; a < P; ++a)
            for (int b = a + 1; b < P; ++b) pr.pairs.push_back(std::make_pair(a, b));
        for (size_t i = 0; i < pr.pairs.size(); ++i)
            pr.pair_names.push_back(as<std::string>(ctx.trait_names[pr.pairs[i].first]) + "_" +
                                    as<std::string>(ctx.trait_names[pr.pairs[i].second]));
        Rcout << "Co-heritability: " << pr.pairs.size() << " trait pair(s); vg = genetic covariance, "
              << "h2 = vg / sqrt(Vp1 Vp2), p_spa two-sided.\n";
    }
    // ---- binary co-heritability inputs (exact conditional saddlepoint) ---------
    if (coher && binary) {
        if (binary_raw.isNull())
            stop("binary = TRUE with coher = TRUE needs binary_raw (0/1 matrix, one column per trait)");
        Rcpp::NumericMatrix braw(binary_raw.get());
        if (braw.ncol() != P) stop("binary_raw must have one column per phenotype");
        const int n = ctx.n_inds;
        ctx.braw.assign(P, std::vector<unsigned char>(n, 0));
        ctx.prev.assign(P, 0.0); ctx.bdelta.assign(P, 0.0);
        for (int t = 0; t < P; ++t) {
            double sumb = 0.0, bty = 0.0;
            for (int i = 0; i < n; ++i) {
                const double v = braw(ctx.pheno_keep[i], t);
                if (!(v == 0.0 || v == 1.0)) stop("binary_raw must contain only 0 and 1 (column %d)", t + 1);
                ctx.braw[t][i] = (unsigned char) v; sumb += v; bty += v * ctx.Y(i, t);
            }
            ctx.prev[t] = sumb / n;
            if (!(bty > 0.0)) stop("binary_raw column %d does not match the phenotype", t + 1);
            ctx.bdelta[t] = (double)(n - 1) / bty;      // cancels in the p-value anyway
        }
        Rcout << "Binary co-heritability: exact conditional Bernoulli saddlepoint; prevalences:";
        for (int t = 0; t < P; ++t) Rcout << " " << ctx.prev[t];
        Rcout << "\n";
    }
    if (binary && !coher)
        Rcout << "Binary phenotype: 4th-cumulant variance correction + six-cumulant saddlepoint.\n";

    if (!Rf_isNull(chr)) {
        Rcpp::CharacterVector cs(Rf_coerceVector(chr, STRSXP));
        if (cs.size() != 1) stop("chr must be a single value");
        pr.chr = Rcpp::as<std::string>(cs[0]);
        bool found = false;
        for (size_t i = 0; i < ctx.chr_order.size() && !found; ++i)
            if (strip_chr(ctx.chr_order[i]) == strip_chr(pr.chr)) found = true;
        if (!found) stop("chr = '" + pr.chr + "' not found in " + filename + ".bim");
        Rcout << "Restricting to chromosome " << pr.chr << "\n";
    }
    return run_all(ctx, pr);
}
