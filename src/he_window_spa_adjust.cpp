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
#include <climits>
#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <unordered_map>
#include <functional>

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
    std::vector<std::string> id;       // SNP name (column 2): used to match KVIK summaries
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
        out.id.push_back(rsid);
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
    std::vector<int> bim;          // per-column .bim index, -1 for a pseudo-marker (collapsing only)
    // the ORIGINAL variants of the range (collapsing only): MAF and category, so that
    // ACAT-V and the alpha-model weights can still see every variant individually
    std::vector<float> omaf;
    std::vector<int> ocat;
    // exact binary SPA: per column, the carriers of the MINOR allele and their minor-allele
    // dosage (1, 2, or more for a collapsed pseudo-marker). Filled only when requested
    // (car_max > 0) and only for columns with at most car_max carriers; empty = dense.
    std::vector< std::vector<int> > car;
    std::vector< std::vector<unsigned char> > dos;
};

// Carriers of one RAW column (codes 0/1/2, -1 = missing). minor = true: the column is a
// collapsed pseudo-marker whose values already ARE minor-allele dosages.
static void capture_carriers(const float* col, int n, int car_max, bool minor,
                             std::vector<int>& idx, std::vector<unsigned char>& d) {
    idx.clear(); d.clear();
    if (car_max <= 0) return;
    bool flip = false;
    if (!minor) {
        double s = 0.0; int nv = 0;
        for (int i = 0; i < n; ++i) if (col[i] >= 0.f) { s += col[i]; ++nv; }
        flip = (nv > 0 && s > (double) nv);              // coded allele is the major one
    }
    int cnt = 0;
    for (int i = 0; i < n; ++i) {
        const float g = col[i];
        if (g < 0.f) continue;                            // missing = non-carrier
        const int md = (int) std::lround(flip ? 2.f - g : g);
        if (md > 0) { if (++cnt > car_max) { idx.clear(); d.clear(); return; }
                      idx.push_back(i); d.push_back((unsigned char) std::min(md, 255)); }
    }
}

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
                          int collapse_mac = 0, int collapse_n = 5, int car_max = 0) {
    cell = Cell(); cell.i0 = -1;
    int m = i1 - i0;
    int n_keep = (int) keep.size();
    if (m <= 0) { cell.X = GenoMat(n_keep, 0); return true; }
    GenoMat block;                                  // n_keep x m, RAW codes
    if (!br.read_raw(keep, i0, i1, block)) return false;

    // -------- no collapsing: original path (columns map 1:1 to the .bim) --------
    if (!(collapse_mac > 0 && cat_id != 0 && collapse_n >= 1)) {
        if (car_max > 0) {
            cell.car.resize(m); cell.dos.resize(m);
            for (int j = 0; j < m; ++j)
                capture_carriers(block.col(j).data(), n_keep, car_max, false, cell.car[j], cell.dos[j]);
        }
        standardize_capture_maf(block, cell.maf);   // standardize in place
        cell.X = std::move(block);
        cell.i0 = i0;
        return true;
    }

    // -------- collapsing (as SAIGE-GENE+): every variant with MAC <= collapse_mac --
    // is added (sum of MINOR-allele dosages) into ONE pseudo-marker per category of
    // the range. Variants with MAC > collapse_mac stay individual. A pseudo-marker has
    // MAF 0 here, so it is exempt from alpha weighting (its frequency is not defined),
    // and .bim index -1. The original variants are kept in omaf / ocat.
    // collapse_n is no longer used (kept in the signature for compatibility).
    (void) collapse_n;
    int ncat_all = (int) cat_id->size();
    std::vector< std::vector<float> > out_cols;   // each length n_keep (raw scale)
    std::vector<int>  out_cat, out_bim;           // per output column: category id, .bim index
    std::map<int, std::vector<float> > grp;       // category -> pseudo-marker being built
    cell.omaf.assign(m, 0.f); cell.ocat.assign(m, -1);

    for (int j = 0; j < m; ++j) {
        // fold to minor-allele dosage and compute MAC (missing excluded from MAC,
        // treated as 0 dosage for the burden)
        double s = 0.0; int nv = 0;
        for (int i = 0; i < n_keep; ++i) { float g = block(i, j); if (g >= 0.f) { s += g; ++nv; } }
        int si = (int) (s + 0.5);
        int mac = (nv > 0) ? std::min(si, 2 * nv - si) : 0;
        bool flip = (nv > 0 && s > (double) nv);          // coded allele is major
        int cid = (i0 + j >= 0 && i0 + j < ncat_all) ? (*cat_id)[i0 + j] : -1;
        double f = (nv > 0) ? s / (2.0 * nv) : 0.0; if (f > 0.5) f = 1.0 - f;
        cell.omaf[j] = (float) f; cell.ocat[j] = cid;

        if (mac > collapse_mac) {
            std::vector<float> col(n_keep);
            for (int i = 0; i < n_keep; ++i) { float g = block(i, j); col[i] = (g < 0.f) ? -1.f : g; }
            out_cols.push_back(std::move(col)); out_cat.push_back(cid); out_bim.push_back(i0 + j);
        } else if (mac > 0) {
            std::vector<float>& acc = grp[cid];
            if (acc.empty()) acc.assign(n_keep, 0.f);
            for (int i = 0; i < n_keep; ++i) {
                float g = block(i, j);
                if (g >= 0.f) acc[i] += (flip ? (2.f - g) : g);   // additive minor dosage
            }
        }
    }
    const int n_ind = (int) out_cols.size();
    for (std::map<int, std::vector<float> >::iterator it = grp.begin(); it != grp.end(); ++it) {
        out_cols.push_back(it->second); out_cat.push_back(it->first); out_bim.push_back(-1);
    }

    int nout = (int) out_cols.size();
    if (car_max > 0) {
        cell.car.resize(nout); cell.dos.resize(nout);
        for (int k = 0; k < nout; ++k)
            capture_carriers(out_cols[k].data(), n_keep, car_max, k >= n_ind, cell.car[k], cell.dos[k]);
    }
    GenoMat X(n_keep, nout);
    for (int k = 0; k < nout; ++k)
        for (int i = 0; i < n_keep; ++i) X(i, k) = out_cols[k][i];
    standardize_capture_maf(X, cell.maf);
    for (int k = n_ind; k < nout; ++k) cell.maf[k] = 0.f;   // pseudo-markers: exempt from alpha
    cell.X = std::move(X);
    cell.cat = out_cat;
    cell.bim = out_bim;
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

// All p-values here are TWO-SIDED: P(|X - mean| >= |s_obs - mean|), reported as
// 2 * (tail in the direction of the observation). Both the univariate variance
// component and the co-heritability use this convention.
static QuadSpaResult quad_spa_solve(
    double s_obs_in, const std::vector<double>& eig_in, double eig_rep_in, double n_rep,
    int max_iter = 100, double tol = 1e-8,
    double lam2_in = 0.0, double nu2 = 0.0,
    double lam3_in = 0.0, double nu3 = 0.0,
    int tail = 0                     // 0 = two-sided (from the median), 1 = upper tail P(X >= s)
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
        res.p = (tail == 1) ? 0.5 * std::erfc(z / std::sqrt(2.0)) : std::erfc(std::abs(z) / std::sqrt(2.0));
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
    const double lead = (w >= 0) ? 0.5 * std::erfc(w / std::sqrt(2.0)) : Phi_w;
    const double corr = phi_w * (1.0 / u - 1.0 / w);
    double p_one = (w >= 0) ? lead + corr : lead - corr;

    // ---- VALIDITY OF THE EXPANSION -- do not replace with a clamp ----------
    // Lugannani-Rice is an ASYMPTOTIC expansion, not a bound. When the fitted
    // CGF is badly conditioned -- which is what the six-cumulant fit produces
    // from the enormous high-order cumulants of a very low prevalence binary
    // trait -- the saddlepoint lands far out, |u| >> |w|, and the correction
    // term swamps the leading term and drives p_one negative.
    //
    // Clamping that to 0 and returning converged = true reported a statistic
    // 1.66 SE from zero as p = 0. Reject instead: the caller then falls back
    // to the 4-cumulant fit, then the eigenvalue path, then the Wald p, all of
    // which are honest. A rejected saddlepoint costs a little power; a silent
    // p = 0 costs the whole calibration.
    if (!std::isfinite(p_one) || p_one <= 0.0 || p_one >= 1.0) return res;
    if (std::abs(corr) > std::abs(lead)) return res;   // 1st-order term exceeds the 0th

    // Removable singularity at w ~ 0 (distinct from the t ~ 0 test above: the
    // expansion is unstable in w even when t is not small).
    if (std::abs(w) < 1e-4) {
        const double z = (s_obs - mean0) / std::sqrt(var0);
        res.p = (tail == 1) ? 0.5 * std::erfc(z / std::sqrt(2.0)) : std::erfc(std::abs(z) / std::sqrt(2.0));
        res.converged = true;
        return res;
    }
    if (tail == 1) {                   // ONE-SIDED (variance components): only large values count
        res.p = (w >= 0) ? p_one : 1.0 - p_one;
        res.converged = true;
        return res;
    }

    // TWO-SIDED, measured from the MEDIAN:  p = 2 min{F(x), 1 - F(x)}.
    // p_one is the tail on the observation's side of the MEAN. For a skewed
    // null (small K) the median lies below the mean, so for observations in
    // between, p_one > 0.5 and the old 2 * p_one was capped at exactly 1 --
    // 13-18% of tests for 1-2 SNP bins, each of which then dragged an ACAT
    // combination to ~1. 1 - p_one is the other tail, so taking the smaller
    // one is exact, uniform under the null, and leaves both tails unchanged.
    res.p = 2.0 * std::min(p_one, 1.0 - p_one);
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
    // Same two rules as quad_spa_solve: an expansion that leaves (0,1) has
    // failed and is rejected (the caller falls back to the Wald p) rather than
    // clamped to an exact 0; and the two-sided p is measured from the median.
    if (!std::isfinite(p_one) || p_one <= 0.0 || p_one >= 1.0) return res;
    res.p = 2.0 * std::min(p_one, 1.0 - p_one);
    res.converged = true;
    return res;
}


// ===========================================================================
// EXACT PER-VARIANT SADDLEPOINT FOR BINARY TRAITS  (spa_used = 5)
// ===========================================================================
// Every estimate is a quadratic form in the per-column scores u_j = v_j' y:
//     Q = sum_j a_j u_j^2 + C ,   a_j = Tinv(c, comp(j)) g_comp(j),  C = Tinv(c, env) y'y
// For a 0/1 trait, y_i = lo_i + delta b_i with b_i in {0,1}. For a SPARSE column (few
// expected case carriers) u_j is, up to a fixed shift, delta * kappa_j * S_j, where S_j is
// the (dosage-weighted) number of cases among the column's carriers -- a Poisson-binomial
// count whose distribution is computed EXACTLY. Sparse columns rarely share carriers, so
// their terms are treated as independent. DENSE columns (common variants, or rare
// variants with many expected case carriers) are close to Gaussian; their joint
// contribution is a non-central chi-square mixture from their exact covariance.
//     K(t) = C t + sum_sparse log E exp(t a_j u_j^2) + sum_dense K_ncchi2(t)
// y'y is treated as fixed (the null conditions on the number of cases).
// The CGF is evaluated stably with log-sum-exp; the saddlepoint uses Lugannani-Rice
// with the same guards as quad_spa_solve. Returns BOTH tails (upper = P(Q >= q)).
// ---------------------------------------------------------------------------
struct ExactTerm { std::vector<double> v, lp; int obs = -1; };   // support values, log probabilities,
                                                                 // index of the observed value (-1 = unknown)

struct ExactSpaResult { double p_up, p_lo; bool converged; bool exact = false; };   // exact: enumeration / atom

// Many small Bernoulli terms, grouped: w independent copies of c * Bernoulli(p) per group.
// Used for the part of a binary co-heritability statistic that is NOT kept as individual
// exact terms; its exact CGF replaces a Gaussian approximation, which is badly anti-
// conservative when p is small (a sum of Bernoulli(0.01) terms is Poisson-like, skewed).
struct BernRest { std::vector<double> c, p, w; };

static void exact_cgf(double t, const std::vector<ExactTerm>& T, const std::vector<double>& lam,
                      const std::vector<double>& nc2, double C,
                      double& K, double& K1, double& K2, const BernRest* br = 0) {
    K = C * t; K1 = C; K2 = 0.0;
    if (br) for (size_t g = 0; g < br->c.size(); ++g) {
        const double c = br->c[g], p = br->p[g], w = br->w[g], x = t * c;
        double lk, qq;                                   // log(1 - p + p e^x), tilted probability
        if (x > 0.0) { const double d = p + (1.0 - p) * std::exp(-x); lk = x + std::log(d); qq = p / d; }
        else         { const double e = p * std::exp(x); lk = std::log1p(e - p); qq = e / (1.0 - p + e); }
        K += w * lk; K1 += w * c * qq; K2 += w * c * c * qq * (1.0 - qq);
    }
    for (size_t j = 0; j < T.size(); ++j) {
        const std::vector<double>& v = T[j].v; const std::vector<double>& lp = T[j].lp;
        double mx = -1e300;
        for (size_t s = 0; s < v.size(); ++s) mx = std::max(mx, lp[s] + t * v[s]);
        double z = 0.0, m1 = 0.0, m2 = 0.0;
        for (size_t s = 0; s < v.size(); ++s) {
            const double e = std::exp(lp[s] + t * v[s] - mx);
            z += e; m1 += e * v[s]; m2 += e * v[s] * v[s];
        }
        m1 /= z; m2 /= z;
        K += mx + std::log(z); K1 += m1; K2 += std::max(0.0, m2 - m1 * m1);
    }
    for (size_t k = 0; k < lam.size(); ++k) {
        const double l = lam[k], d = 1.0 - 2.0 * l * t;
        K  += -0.5 * std::log(d) + l * nc2[k] * t / d;
        K1 += l / d + l * nc2[k] / (d * d);
        K2 += 2.0 * l * l / (d * d) + 4.0 * l * l * nc2[k] / (d * d * d);
    }
}

static ExactSpaResult exact_spa_solve(double q_in, std::vector<ExactTerm> T, std::vector<double> lam,
                                      const std::vector<double>& nc2, double C_in, const BernRest* br_in = 0) {
    ExactSpaResult res; res.p_up = res.p_lo = NA_REAL; res.converged = false;
    // ---- scale normalisation (as quad_spa_solve): work where Var(Q) = 1 --------------
    double mean = C_in, var = 0.0, qmax = C_in, qmin = C_in;
    for (size_t j = 0; j < T.size(); ++j) {
        double m1 = 0.0, m2 = 0.0, vmx = -1e300, vmn = 1e300;
        for (size_t s = 0; s < T[j].v.size(); ++s) {
            const double pr = std::exp(T[j].lp[s]);
            m1 += pr * T[j].v[s]; m2 += pr * T[j].v[s] * T[j].v[s];
            vmx = std::max(vmx, T[j].v[s]); vmn = std::min(vmn, T[j].v[s]);
        }
        mean += m1; var += std::max(0.0, m2 - m1 * m1); qmax += vmx; qmin += vmn;
    }
    bool unb_up = false, unb_lo = false;
    for (size_t k = 0; k < lam.size(); ++k) {
        mean += lam[k] * (1.0 + nc2[k]); var += 2.0 * lam[k] * lam[k] * (1.0 + 2.0 * nc2[k]);
        if (lam[k] > 0) unb_up = true; if (lam[k] < 0) unb_lo = true;
    }
    BernRest brs; const BernRest* br = 0;
    if (br_in && !br_in->c.empty()) {
        brs = *br_in; br = &brs;
        for (size_t g = 0; g < brs.c.size(); ++g) {
            mean += brs.w[g] * brs.c[g] * brs.p[g]; var += brs.w[g] * brs.c[g] * brs.c[g] * brs.p[g] * (1.0 - brs.p[g]);
        }
        unb_up = unb_lo = true;              // the atom shortcut below covers the terms T only
    }
    if (!(var > 0.0) || !std::isfinite(var)) return res;
    const double sc = 1.0 / std::sqrt(var);
    for (size_t j = 0; j < T.size(); ++j) for (size_t s = 0; s < T[j].v.size(); ++s) T[j].v[s] *= sc;
    for (size_t k = 0; k < lam.size(); ++k) lam[k] *= sc;
    if (br) for (size_t g = 0; g < brs.c.size(); ++g) brs.c[g] *= sc;
    const double C = C_in * sc, q = q_in * sc, mu0 = mean * sc;
    qmax *= sc; qmin *= sc;

    // ---- observation at the edge of a bounded support: exact atom probability ---------
    const double tol = 1e-9 * std::max(1.0, std::abs(q));
    if (!unb_up && q >= qmax - tol) {
        double lp = 0.0;
        for (size_t j = 0; j < T.size(); ++j) {
            double vmx = -1e300; for (size_t s = 0; s < T[j].v.size(); ++s) vmx = std::max(vmx, T[j].v[s]);
            double pr = 0.0; for (size_t s = 0; s < T[j].v.size(); ++s) if (T[j].v[s] >= vmx - tol) pr += std::exp(T[j].lp[s]);
            lp += std::log(std::max(pr, 1e-300));
        }
        res.p_up = std::max(std::exp(lp), 1e-300); res.p_lo = 1.0; res.converged = true; res.exact = true; return res;
    }
    if (!unb_lo && q <= qmin + tol) {
        double lp = 0.0;
        for (size_t j = 0; j < T.size(); ++j) {
            double vmn = 1e300; for (size_t s = 0; s < T[j].v.size(); ++s) vmn = std::min(vmn, T[j].v[s]);
            double pr = 0.0; for (size_t s = 0; s < T[j].v.size(); ++s) if (T[j].v[s] <= vmn + tol) pr += std::exp(T[j].lp[s]);
            lp += std::log(std::max(pr, 1e-300));
        }
        res.p_lo = std::max(std::exp(lp), 1e-300); res.p_up = 1.0; res.converged = true; res.exact = true; return res;
    }

    // ---- purely discrete null (no Gaussian part): enumerate the distribution -----------
    // The support of the sum is built term by term, merging equal values. Combinations
    // with probability < 1e-14 are dropped and their total mass is ADDED to the upper-
    // tail p (and to the lower-tail p), so the result is an exact upper bound that is
    // tight for p >> 1e-14. Used when the product of the term supports is <= 200,000;
    // abandoned (-> saddlepoint) if the merged support exceeds 50,000 points.
    double prod_support = 1.0;                 // cheap pre-check of the enumeration size
    for (size_t j = 0; j < T.size() && prod_support <= 2e5; ++j) {
        int cnt = 0; for (size_t s2 = 0; s2 < T[j].lp.size(); ++s2) if (T[j].lp[s2] >= -32.2) ++cnt;   // p >= 1e-14
        prod_support *= std::max(1, cnt);
    }
    if (lam.empty() && !br && prod_support <= 2e5) {
        std::vector< std::pair<double,double> > dist(1, std::make_pair(C, 1.0)), nx;
        bool ok = true; double dropped = 0.0;
        for (size_t j = 0; j < T.size() && ok; ++j) {
            nx.clear(); nx.reserve(dist.size() * T[j].v.size());
            for (size_t a2 = 0; a2 < dist.size(); ++a2)
                for (size_t s2 = 0; s2 < T[j].v.size(); ++s2) {
                    const double pr = dist[a2].second * std::exp(T[j].lp[s2]);
                    if (pr >= 1e-14) nx.push_back(std::make_pair(dist[a2].first + T[j].v[s2], pr));
                    else dropped += pr;
                }
            std::sort(nx.begin(), nx.end());
            dist.clear();
            for (size_t k = 0; k < nx.size(); ++k) {
                if (!dist.empty() && nx[k].first - dist.back().first <= 1e-12) dist.back().second += nx[k].second;
                else dist.push_back(nx[k]);
            }
            if (dist.size() > 50000) ok = false;
        }
        if (ok) {
            double pu = 0.0, pl = 0.0;
            for (size_t k = 0; k < dist.size(); ++k) {
                if (dist[k].first >= q - tol) pu += dist[k].second;
                if (dist[k].first <= q + tol) pl += dist[k].second;
            }
            res.p_up = std::min(1.0, std::max(pu + dropped, 1e-300));
            res.p_lo = std::min(1.0, std::max(pl + dropped, 1e-300));
            res.converged = true; res.exact = true; return res;
        }
    }

    // ---- saddlepoint: K'(t) = q, bracketed Newton (K' is increasing) -----------------
    double lmax = 0.0, lmin = 0.0;
    for (size_t k = 0; k < lam.size(); ++k) { lmax = std::max(lmax, lam[k]); lmin = std::min(lmin, lam[k]); }
    double hi = (lmax > 1e-14) ? (1.0 - 1e-9) / (2.0 * lmax) :  1e4;
    double lo = (lmin < -1e-14) ? (1.0 - 1e-9) / (2.0 * lmin) : -1e4;
    double K, K1, K2, t = 0.0;
    exact_cgf(0.0, T, lam, nc2, C, K, K1, K2, br);
    double a = lo, b = hi;
    bool conv = false;
    for (int it = 0; it < 300; ++it) {
        exact_cgf(t, T, lam, nc2, C, K, K1, K2, br);
        const double f = K1 - q;
        if (std::abs(f) < 1e-10 * std::max(1.0, std::abs(q - mu0))) { conv = true; break; }
        if (f > 0) b = t; else a = t;
        double tn = (K2 > 1e-300) ? t - f / K2 : 0.5 * (a + b);
        if (!(tn > a && tn < b) || !std::isfinite(tn)) tn = 0.5 * (a + b);
        if (std::abs(b - a) < 1e-14 * std::max(1.0, std::abs(t))) { conv = true; t = tn; break; }
        t = tn;
    }
    if (!conv) return res;
    exact_cgf(t, T, lam, nc2, C, K, K1, K2, br);
    if (!(K2 > 0.0) || !std::isfinite(K2)) return res;

    double up, dn;
    if (std::abs(t) < 1e-7) {
        const double z = (q - mu0);                        // Var = 1 after scaling
        up = 0.5 * std::erfc(z / std::sqrt(2.0)); dn = 1.0 - up;
    } else {
        const double w = ((t > 0) ? 1.0 : -1.0) * std::sqrt(std::max(0.0, 2.0 * (t * q - K)));
        const double u = t * std::sqrt(K2);
        if (std::abs(w) < 1e-4) {
            const double z = (q - mu0);
            up = 0.5 * std::erfc(z / std::sqrt(2.0)); dn = 1.0 - up;
        } else {
            const double phi = std::exp(-0.5 * w * w) / std::sqrt(2.0 * M_PI);
            const double corr = phi * (1.0 / u - 1.0 / w);
            if (w >= 0) {
                const double lead = 0.5 * std::erfc(w / std::sqrt(2.0));
                if (std::abs(corr) > std::abs(lead)) return res;
                up = lead + corr; dn = 1.0 - up;
            } else {
                // Below the mean. Near the lower bound of a discrete null (few sparse terms,
                // no Gaussian part) the expansion can break down; there the UPPER-tail p is
                // > 0.5 anyway, so the normal approximation is used instead of failing.
                const double lead = 0.5 * std::erfc(-w / std::sqrt(2.0));
                if (std::abs(corr) > std::abs(lead)) {
                    dn = 0.5 * std::erfc(-(q - mu0) / std::sqrt(2.0)); up = 1.0 - dn;
                } else {
                    dn = lead - corr; up = 1.0 - dn;
                }
            }
        }
    }
    if (!std::isfinite(up) || !std::isfinite(dn) || up <= 0.0 || dn <= 0.0 || up >= 1.0 || dn >= 1.0) return res;
    res.p_up = up; res.p_lo = dn; res.converged = true;
    return res;
}

// ---------------------------------------------------------------------------
// Lattice convolution of the same null (used when the saddlepoint fails or p is small).
// Why: a sparse term a_j u_j^2 is the SQUARE of a (Poisson-)binomial count, so its
// moment generating function grows super-exponentially -- it is finite only because the
// support is bounded. For a category made of variants with tens of carriers each (no
// singletons; e.g. after collapsing) the saddlepoint then lands among atoms with
// probability < 1e-30 and Lugannani-Rice is meaningless (it fails, or is far off).
// Here the sparse part D = sum_j (v_j - min_j) is computed on a lattice of M bins
// (each atom rounded to the nearest bin; errors average out over terms), with an
// absorbing overflow bin at the top. The dense part G (weighted non-central chi-square)
// is continuous and its tail is accurate by saddlepoint, so
//     P(Q >= q) = sum_k P(D = k h) P(G >= q - C - sum_j min_j - k h).
// Without a dense part the threshold is the observed lattice point (sum of the rounded
// observed atoms), so the observed configuration itself is always counted, as in the
// enumeration. Atoms with probability < 1e-18 are dropped and their mass is ADDED to
// both tails (an upper bound).
// ---------------------------------------------------------------------------
static ExactSpaResult exact_grid_solve(double q, const std::vector<ExactTerm>& T,
                                       const std::vector<double>& lam, const std::vector<double>& nc2,
                                       double C, int M = 8192, const BernRest* br = 0) {
    ExactSpaResult res; res.p_up = res.p_lo = NA_REAL; res.converged = false;
    const double LPMIN = std::log(1e-18);
    if (br && br->c.empty()) br = 0;
    const bool dense = !lam.empty() || br;
    // per term: kept atoms, minimum, observed value
    double sum_mn = 0.0, Dobs = 0.0, dropped = 0.0; bool have_obs = !dense;
    std::vector<double> mn(T.size(), 0.0);
    for (size_t j = 0; j < T.size(); ++j) {
        double m = 1e300;
        for (size_t s = 0; s < T[j].v.size(); ++s)
            if (T[j].lp[s] >= LPMIN || (int) s == T[j].obs) m = std::min(m, T[j].v[s]);
        if (!(m < 1e300)) return res;
        mn[j] = m; sum_mn += m;
        if (T[j].obs < 0 || T[j].obs >= (int) T[j].v.size()) have_obs = false;
        else Dobs += T[j].v[T[j].obs] - m;
    }
    const double Dq = q - C - sum_mn;                    // threshold for D (+ G)
    // dense part: moments, support sign
    double mG = 0.0, vG = 0.0; bool allpos = true, allneg = true;
    for (size_t k = 0; k < lam.size(); ++k) {
        mG += lam[k] * (1.0 + nc2[k]); vG += 2.0 * lam[k] * lam[k] * (1.0 + 2.0 * nc2[k]);
        if (lam[k] < 0) allpos = false;
        if (lam[k] > 0) allneg = false;
    }
    double Gmin = 0.0, Gmax = 0.0;                       // support of the Bernoulli part
    if (br) for (size_t g = 0; g < br->c.size(); ++g) {
        const double c = br->c[g], p = br->p[g], w = br->w[g];
        mG += w * c * p; vG += w * c * c * p * (1.0 - p);
        if (c < 0) { allpos = false; Gmin += w * c; } else Gmax += w * c;
        if (c > 0) allneg = false;
    }
    const bool bounded = lam.empty();                    // G = Bernoulli part only
    const double sdG = std::sqrt(std::max(vG, 0.0));
    double R;
    if (!dense) R = have_obs ? Dobs : Dq;
    else        R = Dq - (allpos ? std::max(0.0, mG - 12.0 * sdG) : (mG - 12.0 * sdG));
    if (!dense && !(R > 0.0)) {                          // observation at the minimum of D
        res.p_up = 1.0; res.p_lo = 1.0; res.converged = true; return res;
    }
    double scale = 0.0;                                  // a natural size, for degenerate R
    for (size_t j = 0; j < T.size(); ++j) for (size_t s = 0; s < T[j].v.size(); ++s)
        if (T[j].lp[s] >= LPMIN) scale = std::max(scale, T[j].v[s] - mn[j]);
    if (!(R > 0.0)) R = 1e-9 * std::max(scale, 1e-300);
    const double h = R / M;
    long Kobs = M;
    if (!dense && have_obs) {
        Kobs = 0;
        for (size_t j = 0; j < T.size(); ++j) Kobs += std::llround((T[j].v[T[j].obs] - mn[j]) / h);
    }
    const long Mc = std::max<long>(M, Kobs) + 1;         // bins 0..Mc-1, then overflow
    std::vector<double> pm(Mc, 0.0), nx(Mc, 0.0);
    pm[0] = 1.0; long top = 0; double ovf = 0.0;
    std::vector<long> ak; std::vector<double> ap;
    for (size_t j = 0; j < T.size(); ++j) {
        ak.clear(); ap.clear();
        for (size_t s = 0; s < T[j].v.size(); ++s) {
            const double pr = std::exp(T[j].lp[s]);
            if (T[j].lp[s] < LPMIN && (int) s != T[j].obs) { dropped += pr; continue; }
            ak.push_back(std::llround((T[j].v[s] - mn[j]) / h)); ap.push_back(pr);
        }
        long ntop = 0;
        std::fill(nx.begin(), nx.end(), 0.0);
        for (long k = 0; k <= top; ++k) {
            const double pk = pm[k]; if (pk == 0.0) continue;
            for (size_t a = 0; a < ak.size(); ++a) {
                const long idx = k + ak[a];
                if (idx >= Mc) ovf += pk * ap[a];
                else { nx[idx] += pk * ap[a]; if (idx > ntop) ntop = idx; }
            }
        }
        pm.swap(nx); top = ntop;
    }
    double pu = ovf, pl = 0.0;
    if (!dense) {
        for (long k = 0; k <= top; ++k) { if (k >= Kobs) pu += pm[k]; if (k <= Kobs) pl += pm[k]; }
    } else {
        // tail of G at x: saddlepoint of the non-central chi-square mixture; cached on a
        // coarse stride and interpolated on the log scale (G's tail is smooth)
        const std::vector<ExactTerm> none;
        auto pg = [&](double x, double& up, double& lo) {
            if (allpos && x <= 0.0) { up = 1.0; lo = 0.0; return; }     // G >= 0
            if (allneg && x >= 0.0) { up = 0.0; lo = 1.0; return; }     // G <= 0
            if (bounded && x > Gmax) { up = 0.0; lo = 1.0; return; }
            if (bounded && x < Gmin) { up = 1.0; lo = 0.0; return; }
            ExactSpaResult r = exact_spa_solve(x, none, lam, nc2, 0.0, br);
            if (r.converged) { up = r.p_up; lo = r.p_lo; }
            else { const double z = (x - mG) / std::max(sdG, 1e-300);
                   up = 0.5 * std::erfc(z / std::sqrt(2.0)); lo = 1.0 - up; }
        };
        const long stride = 8;
        double lu0 = 0, ll0 = 0, lu1 = 0, ll1 = 0; long k0 = -1;
        for (long k = 0; k <= top; ++k) {
            const double pk = pm[k]; if (pk < 1e-300) continue;
            const long kb = (k / stride) * stride;
            if (kb != k0) {
                double u0, l0, u1, l1;
                pg(Dq - kb * h, u0, l0); pg(Dq - (kb + stride) * h, u1, l1);
                lu0 = std::log(std::max(u0, 1e-300)); ll0 = std::log(std::max(l0, 1e-300));
                lu1 = std::log(std::max(u1, 1e-300)); ll1 = std::log(std::max(l1, 1e-300));
                k0 = kb;
            }
            const double f = (double) (k - kb) / stride;
            pu += pk * std::exp(lu0 + f * (lu1 - lu0));
            pl += pk * std::exp(ll0 + f * (ll1 - ll0));
        }
        { double u, l; pg(Dq - Mc * h, u, l); pl += ovf * l; }
    }
    res.p_up = std::min(1.0, std::max(pu + dropped, 1e-300));
    res.p_lo = std::min(1.0, std::max(pl + dropped, 1e-300));
    res.converged = std::isfinite(res.p_up) && std::isfinite(res.p_lo);
    return res;
}

// ---------------------------------------------------------------------------
// Lattice convolution for BINARY CO-HERITABILITY: Q = C + sum_j v_j b_j + G, with b_j ~
// Bernoulli(p_j) two-point terms (the carriers) and G the Bernoulli rest (BernRest, tail by
// saddlepoint). Unlike exact_grid_solve (which shifts every term to its minimum, fine for
// the squared counts of single traits), each term is centred at its MOST LIKELY value, so
// that the lattice only has to resolve the few terms that are "on": with hundreds of terms
// with the same negative coefficient, shifting to the minimum puts the bulk of the
// distribution far from zero, and the rounding errors of all terms then add up instead of
// cancelling. The lattice spans [-R-, R+], where R+ (R-) is the sum of the k largest positive
// (negative) increments and k is chosen such that P(more than k positive (negative) terms
// are on) < 1e-16 (from the exact distribution of the counts); that mass is added to both p.
// ---------------------------------------------------------------------------
static ExactSpaResult coher_grid_solve(double q, const std::vector<double>& v, const std::vector<double>& pr,
                                       const BernRest* br, double C, int M = 8192) {
    ExactSpaResult res; res.p_up = res.p_lo = NA_REAL; res.converged = false;
    if (br && br->c.empty()) br = 0;
    double base = C;
    std::vector<double> dp, pp, dn, pn;             // positive / negative increments, probabilities
    for (size_t j = 0; j < v.size(); ++j) {
        double d = v[j], p = pr[j];
        if (p > 0.5) { base += d; d = -d; p = 1.0 - p; }
        if (d == 0.0 || p <= 0.0) continue;
        if (d > 0) { dp.push_back(d); pp.push_back(p); } else { dn.push_back(-d); pn.push_back(p); }
    }
    double dropped = 0.0;
    // range of one side: sum of the k largest increments, k from the count distribution
    auto side_range = [&](const std::vector<double>& d, const std::vector<double>& p) -> double {
        if (d.empty()) return 0.0;
        const int J = (int) d.size(), KC = std::min(J, 400);
        std::vector<double> cnt(KC + 2, 0.0); cnt[0] = 1.0;      // P(N = k), k <= KC, last = overflow
        for (int j = 0; j < J; ++j) {
            for (int k = KC + 1; k >= 1; --k) {
                if (k == KC + 1) cnt[k] += cnt[k - 1] * p[j];         // absorbing
                else cnt[k] = cnt[k] * (1.0 - p[j]) + cnt[k - 1] * p[j];
            }
            cnt[0] *= (1.0 - p[j]);
        }
        double tail = cnt[KC + 1]; int k = KC;
        while (k > 0 && tail + cnt[k] < 1e-16) { tail += cnt[k]; --k; }
        dropped += tail;
        std::vector<double> ds(d); std::sort(ds.begin(), ds.end(), std::greater<double>());
        double r = 0.0; for (int i = 0; i < std::min(k, J); ++i) r += ds[i];
        return r;
    };
    const double Rp = side_range(dp, pp), Rn = side_range(dn, pn);
    // the Bernoulli rest G: moments and support
    double mG = 0.0, vG = 0.0, Gmin = 0.0, Gmax = 0.0;
    if (br) for (size_t g = 0; g < br->c.size(); ++g) {
        const double c = br->c[g], p = br->p[g], w = br->w[g];
        mG += w * c * p; vG += w * c * c * p * (1.0 - p);
        if (c < 0) Gmin += w * c; else Gmax += w * c;
    }
    const double x0 = q - base;                     // threshold for D + G
    const double span = Rp + Rn;
    const double sdG = std::sqrt(std::max(vG, 0.0));
    if (!(span > 0.0)) {                            // no discrete part: G only
        if (!br) return res;
        ExactSpaResult r = exact_spa_solve(x0, std::vector<ExactTerm>(), std::vector<double>(), std::vector<double>(), 0.0, br);
        return r;
    }
    const double h = span / M;
    const long L0 = std::lround(Rn / h);            // index of D = 0
    const long NB = M + 1;
    std::vector<double> pm(NB, 0.0), nx(NB, 0.0);
    pm[L0] = 1.0; long lo = L0, hi = L0;
    auto add = [&](double d, double p) {
        const long sh = std::lround(d / h);
        if (sh == 0) return;                        // below lattice resolution: no effect
        std::fill(nx.begin() + std::max(0L, std::min(lo, lo + sh)), nx.begin() + std::min(NB, std::max(hi, hi + sh) + 1), 0.0);
        long nlo = NB, nhi = -1;
        for (long k = lo; k <= hi; ++k) {
            const double pk = pm[k]; if (pk == 0.0) continue;
            nx[k] += pk * (1.0 - p);
            const long t = k + sh;
            if (t < 0 || t >= NB) dropped += pk * p;
            else { nx[t] += pk * p; nlo = std::min(nlo, t); nhi = std::max(nhi, t); }
            nlo = std::min(nlo, k); nhi = std::max(nhi, k);
        }
        for (long k = std::max(0L, std::min(lo, lo + sh)); k < std::min(NB, std::max(hi, hi + sh) + 1); ++k) { pm[k] = nx[k]; }
        lo = nlo; hi = nhi;
    };
    for (size_t j = 0; j < dp.size(); ++j) add(dp[j], pp[j]);
    for (size_t j = 0; j < dn.size(); ++j) add(-dn[j], pn[j]);
    double pu = 0.0, pl = 0.0;
    if (!br) {
        const double tol = 2.0 * h;                 // counts the observed configuration (rounding)
        for (long k = lo; k <= hi; ++k) {
            const double D = (k - L0) * h;
            if (D >= x0 - tol) pu += pm[k];
            if (D <= x0 + tol) pl += pm[k];
        }
    } else {
        const std::vector<ExactTerm> none; const std::vector<double> nolam;
        auto pg = [&](double x, double& up, double& lw) {
            if (x > Gmax) { up = 0.0; lw = 1.0; return; }
            if (x < Gmin) { up = 1.0; lw = 0.0; return; }
            ExactSpaResult r = exact_spa_solve(x, none, nolam, nolam, 0.0, br);
            if (r.converged) { up = r.p_up; lw = r.p_lo; }
            else { const double z = (x - mG) / std::max(sdG, 1e-300);
                   up = 0.5 * std::erfc(z / std::sqrt(2.0)); lw = 1.0 - up; }
        };
        const long stride = 8;
        double lu0 = 0, ll0 = 0, lu1 = 0, ll1 = 0; long k0 = LONG_MIN;
        for (long k = lo; k <= hi; ++k) {
            const double pk = pm[k]; if (pk < 1e-300) continue;
            const long kb = lo + ((k - lo) / stride) * stride;
            if (kb != k0) {
                double u0, l0, u1, l1;
                pg(x0 - (kb - L0) * h, u0, l0); pg(x0 - (kb + stride - L0) * h, u1, l1);
                lu0 = std::log(std::max(u0, 1e-300)); ll0 = std::log(std::max(l0, 1e-300));
                lu1 = std::log(std::max(u1, 1e-300)); ll1 = std::log(std::max(l1, 1e-300));
                k0 = kb;
            }
            const double f = (double) (k - kb) / stride;
            pu += pk * std::exp(lu0 + f * (lu1 - lu0));
            pl += pk * std::exp(ll0 + f * (ll1 - ll0));
        }
    }
    res.p_up = std::min(1.0, std::max(pu + dropped, 1e-300));
    res.p_lo = std::min(1.0, std::max(pl + dropped, 1e-300));
    res.converged = std::isfinite(res.p_up) && std::isfinite(res.p_lo);
    return res;
}

// Distribution of S = sum_i d_i B_i, B_i ~ Bernoulli(pi_i) independent, d_i >= 1 integer.
// Binomial closed form when all d_i = 1 and the pi_i are equal; otherwise a DP.
// Returns log P(S = s), s = 0..smax, trimmed where the probability is negligible.
// Test hook (never set in normal use): capture the exact null model of the first
// exact test of category 0, so a test harness can evaluate the tail function directly.
struct ExactCapture { bool done = false; std::vector<ExactTerm> T; std::vector<double> lam, nc2; double C = 0, q = 0; };
static ExactCapture* g_exact_capture = 0;

static void pb_pmf(const std::vector<unsigned char>& d, const std::vector<double>& pi,
                   std::vector<double>& lp) {
    const int c = (int) d.size();
    bool simple = true; double pmin = 1.0, pmax = 0.0;
    for (int i = 0; i < c; ++i) { if (d[i] != 1) simple = false; pmin = std::min(pmin, pi[i]); pmax = std::max(pmax, pi[i]); }
    if (simple && pmax - pmin <= 1e-9 * std::max(pmax, 1e-12)) {
        const double p = std::min(std::max(0.5 * (pmin + pmax), 1e-12), 1.0 - 1e-12);
        lp.assign(c + 1, 0.0);
        const double l1 = std::log(p), l0 = std::log1p(-p);
        for (int k = 0; k <= c; ++k)
            lp[k] = std::lgamma(c + 1.0) - std::lgamma(k + 1.0) - std::lgamma(c - k + 1.0) + k * l1 + (c - k) * l0;
    } else {
        std::vector<double> pm(1, 1.0), nx;
        for (int i = 0; i < c; ++i) {
            const int di = d[i]; const double p = std::min(std::max(pi[i], 1e-12), 1.0 - 1e-12);
            nx.assign(pm.size() + di, 0.0);
            for (size_t k = 0; k < pm.size(); ++k) { nx[k] += (1.0 - p) * pm[k]; nx[k + di] += p * pm[k]; }
            size_t L = nx.size(); while (L > 1 && nx[L - 1] < 1e-300) --L;
            nx.resize(L); pm.swap(nx);
        }
        lp.resize(pm.size());
        for (size_t k = 0; k < pm.size(); ++k) lp[k] = (pm[k] > 0.0) ? std::log(pm[k]) : -1e300;
    }
    // trim the far upper tail (relative probability < 1e-60): it cannot affect the result
    double mx = -1e300; for (size_t k = 0; k < lp.size(); ++k) mx = std::max(mx, lp[k]);
    size_t L = lp.size(); while (L > 1 && lp[L - 1] < mx - 138.0) --L;
    lp.resize(L);
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
// Regularized incomplete beta I_x(a, b) (continued fraction, Numerical Recipes 6.4), and the
// two-sided Student t p value 2 P(T_nu >= |t|) = I_{nu / (nu + t^2)}(nu / 2, 1 / 2). Accurate in
// the far tail (no cancellation: the direct branch is used when x is small).
static double betacf_(double a, double b, double x) {
    const int MAXIT = 300; const double EPS = 1e-15, FPMIN = 1e-300;
    double qab = a + b, qap = a + 1.0, qam = a - 1.0, c = 1.0, d = 1.0 - qab * x / qap;
    if (std::abs(d) < FPMIN) d = FPMIN; d = 1.0 / d; double h = d;
    for (int m = 1; m <= MAXIT; ++m) {
        const int m2 = 2 * m;
        double aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d; if (std::abs(d) < FPMIN) d = FPMIN;
        c = 1.0 + aa / c;  if (std::abs(c) < FPMIN) c = FPMIN;
        d = 1.0 / d; h *= d * c;
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d; if (std::abs(d) < FPMIN) d = FPMIN;
        c = 1.0 + aa / c;  if (std::abs(c) < FPMIN) c = FPMIN;
        d = 1.0 / d; const double del = d * c; h *= del;
        if (std::abs(del - 1.0) < EPS) break;
    }
    return h;
}
static double ibeta_(double a, double b, double x) {
    if (x <= 0.0) return 0.0; if (x >= 1.0) return 1.0;
    const double lbt = std::lgamma(a + b) - std::lgamma(a) - std::lgamma(b) + a * std::log(x) + b * std::log1p(-x);
    if (x < (a + 1.0) / (a + b + 2.0)) return std::exp(lbt) * betacf_(a, b, x) / a;
    return 1.0 - std::exp(lbt) * betacf_(b, a, 1.0 - x) / b;
}
static double t_two_sided_p(double t, double nu) {
    if (!std::isfinite(nu) || nu > 1e6) return std::erfc(std::abs(t) / std::sqrt(2.0));
    return std::min(1.0, ibeta_(0.5 * nu, 0.5, nu / (nu + t * t)));
}

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
    Eigen::MatrixXd T1, T2;            // tested-component parts of S1 / S2 (L' D L, no e I)
    double t1, t2;                     // the tested-component estimates of both traits
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
    {   // tested-component parts, so that coher_spectrum can shrink ONLY them
        std::vector<double> dt1(C, 0.0), dt2(C, 0.0);
        dt1[c] = d1[c]; dt2[c] = d2[c];
        coher_block(L, W, use_wcache, off, C, K, dt1, 0.0, F.T1);
        coher_block(L, W, use_wcache, off, C, K, dt2, 0.0, F.T2);
        F.t1 = sg1[c]; F.t2 = sg2[c];
    }
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
    // Plug-in estimates need not give a positive definite Omega_m. The usual reason is a
    // NEGATIVE estimate of the tested component for one or both traits (it is kept
    // unconstrained, see coher_fit), which makes S1 or S2 nearly singular in small genes.
    // Repairs, in this order:
    //   1. shrink the NEGATIVE tested-component parts of S1 / S2 toward zero;
    //   2. shrink the GENETIC part of S12 (the background components) toward zero.
    // The scalar parts e1, e2 and e12 are NEVER changed: they also give the n - K repeated
    // eigenvalues below, and changing them in Omega_m only (as an earlier version did when it
    // shrank all of S12) breaks the cancellation between the two, which shifted the mean of
    // the null distribution by up to several standard deviations (p values down to 1e-300
    // for phenotypes with environmental correlation).
    const Eigen::MatrixXd I = Eigen::MatrixXd::Identity(K, K);
    const Eigen::MatrixXd G12 = F.S12 - F.e12 * I;           // genetic part of S12
    double tau = 1.0, shrink = 1.0; bool have = false;
    Eigen::LLT<Eigen::MatrixXd> llt;
    for (int it = 0; it < 40; ++it) {
        Om.topLeftCorner(K, K)     = F.S1 - (F.t1 < 0.0 ? (1.0 - tau) : 0.0) * F.T1;
        Om.bottomRightCorner(K, K) = F.S2 - (F.t2 < 0.0 ? (1.0 - tau) : 0.0) * F.T2;
        Om.topRightCorner(K, K)    = shrink * G12 + F.e12 * I;
        Om.bottomLeftCorner(K, K)  = Om.topRightCorner(K, K);
        llt.compute(Om);
        if (llt.info() == Eigen::Success) { have = true; break; }
        if ((F.t1 < 0.0 || F.t2 < 0.0) && tau > 1e-6) tau *= 0.5;
        else if (shrink > 1e-6) shrink *= 0.5;
        else break;
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
    // alpha values. alphas[0] is the PRIMARY one: it weights the genotypes as
    // read, and every vg / h2 / window estimate is at that value. Any further
    // values are used ONLY for association p-values, which are ACAT-combined
    // across all alphas (see Worker). One value = the original behaviour.
    std::vector<double> alphas;
    // gene definitions (optional). Genes may overlap and nest freely -- they are
    // tested marginally, so the only cost of an overlap is re-reading its SNPs.
    std::vector<std::string> gene_name, gene_chr;
    std::vector<long> gene_start, gene_end;
    // LDAK-KVIK step-2 single-variant results (optional), one per trait, indexed
    // by .bim position: z oriented to the MINOR allele (NaN = not in the file)
    // and the minor allele count 2 n MAF. Drive the ACAT-V component.
    std::vector< std::vector<float> > kz, kmac;
    // collapsing (SAIGE-GENE+ style): variants with MAC <= collapse_mac are summed into one
    // pseudo-marker per category of a tested unit (0 = off). Flanks are never collapsed.
    int collapse_mac = 0;
    // exact per-variant binary saddlepoint (needs the raw 0/1 phenotypes in braw)
    bool bin_exact = false;
    int  car_max = 0;                            // carriers recorded for columns with <= car_max
    double ec_max = 10.0;                        // a column is "sparse" if E[# case carriers] <= ec_max
    std::vector< std::vector<float> > bpi;       // [trait][ind] null P(case) (covariate fit)
    std::vector<double> ncase;                   // observed number of cases per trait
};

// ---------------------------------------------------------------------------
// Cauchy (ACAT) combination. Valid under arbitrary dependence between the
// p-values, which is what we need: the categories and blocks of one gene are
// measured on the same individuals.
//
//   T = sum_j w_j tan((0.5 - p_j) pi),   p = 0.5 - atan(T / sum w_j) / pi
//
// BOTH directions are computed through identities that avoid a cancelling
// subtraction, because gene-based p-values routinely reach 1e-50 and beyond:
//
//   tan((0.5 - p) pi) = cot(pi p) = 1 / tan(pi p)
//       0.5 - p rounds to exactly 0.5 once p < ~1e-17, which would floor every
//       smaller p-value at the same 1.9e-17 -- silently, and only in the tail.
//   0.5 - atan(t)/pi = atan(1/t)/pi   for t > 0
//       the same cancellation on the way back out.
//
// Same class of bug as 1 - Phi(w) in the saddlepoint solvers.
// ---------------------------------------------------------------------------
static double acat_combine(const std::vector<double>& p, const std::vector<double>& w) {
    const double PI = 3.14159265358979323846;
    double T = 0.0, wsum = 0.0;
    for (size_t i = 0; i < p.size(); ++i) {
        if (std::isnan(p[i]) || !(w[i] > 0.0)) continue;
        double pi = p[i];
        if (pi < 1e-300)        pi = 1e-300;        // 1/tan() would overflow
        // Cap at 0.999, not 1 - 1e-15: a p of ~1 contributes ~ -3e14 and wipes
        // out any real signal (p = 1e-5 alongside p = 1 combined to ~1). At
        // 0.999 the contribution is ~ -318, so genuine signal survives, and
        // calibration is unchanged (checked for d = 2, 3, 6). A cap of 1 - 1/d
        // was also tried and inflates P(p < .05) to ~0.054.
        if (pi > 0.999)         pi = 0.999;
        T += w[i] / std::tan(PI * pi);
        wsum += w[i];
    }
    if (!(wsum > 0.0)) return NA_REAL;
    const double t = T / wsum;
    double out = (t > 1.0) ? (std::atan(1.0 / t) / PI) : (0.5 - std::atan(t) / PI);
    if (!(out > 0.0)) out = 1e-300;
    if (out > 1.0)    out = 1.0;
    return out;
}

// ---------------------------------------------------------------------------
// |z| such that the TWO-SIDED normal p-value of z equals p:  erfc(|z|/sqrt2) = p.
// Acklam's rational approximation for the lower-tail quantile, then two Newton
// steps in log space (accurate to ~1e-12 relative, down to p ~ 1e-300).
// Self-contained and thread-safe (no R API), so it can run inside the workers.
// Used to turn a saddlepoint p-value back into a standard error:
//     se_vg = |vg| / z(p)
// ---------------------------------------------------------------------------
static double z_from_p2(double p) {
    if (std::isnan(p)) return NA_REAL;
    if (p >= 1.0) return 0.0;
    if (p < 1e-300) p = 1e-300;
    const double q = 0.5 * p;                                  // upper tail
    static const double a[6] = {-3.969683028665376e+01,  2.209460984245205e+02,
                                -2.759285104469687e+02,  1.383577518672690e+02,
                                -3.066479806614716e+01,  2.506628277459239e+00};
    static const double b[5] = {-5.447609879822406e+01,  1.615858368580409e+02,
                                -1.556989798598866e+02,  6.680131188771972e+01,
                                -1.328068155288572e+01};
    static const double c[6] = {-7.784894002430293e-03, -3.223964580411365e-01,
                                -2.400758277161838e+00, -2.549732539343734e+00,
                                 4.374664141464968e+00,  2.938163982698783e+00};
    static const double d[4] = { 7.784695709041462e-03,  3.224671290700398e-01,
                                 2.445134137142996e+00,  3.754408661907416e+00};
    double x;                                                  // Phi^{-1}(q), q <= 0.5
    if (q < 0.02425) {
        const double r = std::sqrt(-2.0 * std::log(q));
        x = (((((c[0]*r + c[1])*r + c[2])*r + c[3])*r + c[4])*r + c[5]) /
            ((((d[0]*r + d[1])*r + d[2])*r + d[3])*r + 1.0);
    } else {
        const double r = q - 0.5, s2 = r * r;
        x = (((((a[0]*s2 + a[1])*s2 + a[2])*s2 + a[3])*s2 + a[4])*s2 + a[5]) * r /
            (((((b[0]*s2 + b[1])*s2 + b[2])*s2 + b[3])*s2 + b[4])*s2 + 1.0);
    }
    double z = -x;                                             // upper-tail quantile
    for (int it = 0; it < 2; ++it) {                           // Newton on log Q(z) = log q
        const double Q = 0.5 * std::erfc(z / std::sqrt(2.0));
        if (!(Q > 0.0)) break;
        const double phi = std::exp(-0.5 * z * z) / std::sqrt(2.0 * M_PI);
        z += (std::log(Q) - std::log(q)) * Q / phi;
    }
    return z > 0.0 ? z : 0.0;
}

// Standard error implied by an estimate and its two-sided p-value.
static double se_from_p(double est, double p) {
    if (std::isnan(est) || std::isnan(p)) return NA_REAL;
    const double z = z_from_p2(p);
    if (!(z > 1e-8)) return NA_REAL;                           // p ~ 1: SE undefined
    return std::abs(est) / z;
}

// Trim first, THEN strip a chr/Chr/CHR prefix. The trim matters: R's
// as.matrix() on a data.frame with any character column formats the numeric
// columns with format(), which right-aligns them to a common width -- so a
// chromosome column arrives as " 1", " 2", ..., "10" and only the widest label
// would ever match the .bim. Silent, and it looks like a build mismatch.
static std::string strip_chr(const std::string& s) {
    const size_t a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return std::string();
    const size_t b = s.find_last_not_of(" \t\r\n");
    const std::string t = s.substr(a, b - a + 1);
    if (t.size() > 3 && (t.compare(0, 3, "Chr") == 0 || t.compare(0, 3, "chr") == 0 ||
                         t.compare(0, 3, "CHR") == 0)) return t.substr(3);
    return t;
}

// ---------------------------------------------------------------------------
// LDAK-KVIK step-2 single-variant results: the .summaries file
//     Predictor  A1  A2  Z  n  A1Freq
// Z is KVIK's test statistic AFTER its saddlepoint correction (empirical SPA for
// quantitative traits, logistic SPA for binary), i.e. its p-value is exactly
// erfc(|Z|/sqrt2). Variants are matched to the .bim by name. Z is re-oriented
// to the MINOR allele (flipped when A1Freq > 0.5), so that a burden of
// ultra-rare variants adds up effects of the same allele. Main thread only.
// ---------------------------------------------------------------------------
static void kvik_load(const std::string& path, ChunkContext& ctx, int t,
                      const std::unordered_map<std::string, int>& idx) {
    std::ifstream f(path.c_str());
    if (!f.is_open()) stop("Cannot open KVIK summaries file: " + path);
    std::string line;
    if (!std::getline(f, line)) stop("Empty KVIK summaries file: " + path);
    std::istringstream hs(line); std::vector<std::string> hdr; std::string tok;
    while (hs >> tok) hdr.push_back(tok);
    int cP = -1, cZ = -1, cN = -1, cF = -1;
    for (size_t i = 0; i < hdr.size(); ++i) {
        if (hdr[i] == "Predictor") cP = (int) i;
        else if (hdr[i] == "Z")      cZ = (int) i;
        else if (hdr[i] == "n")      cN = (int) i;
        else if (hdr[i] == "A1Freq") cF = (int) i;
    }
    if (cP < 0 || cZ < 0 || cN < 0 || cF < 0)
        stop("KVIK file " + path + " needs columns Predictor, Z, n and A1Freq "
             "(use the .summaries file written by --kvik-step2)");
    const float nanf = std::numeric_limits<float>::quiet_NaN();
    ctx.kz[t].assign(ctx.n_snps, nanf);
    ctx.kmac[t].assign(ctx.n_snps, nanf);
    std::vector<std::string> v(hdr.size());
    long nread = 0, nmatch = 0, nbad = 0;
    auto num = [](const std::string& x, double& out) -> bool {
        char* e = 0; out = std::strtod(x.c_str(), &e);
        return e != x.c_str() && *e == '\0' && std::isfinite(out);
    };
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        std::istringstream ls(line);
        size_t k = 0;
        while (k < v.size() && (ls >> v[k])) ++k;
        if (k < v.size()) { ++nbad; continue; }
        ++nread;
        std::unordered_map<std::string, int>::const_iterator it = idx.find(v[cP]);
        if (it == idx.end()) continue;
        double z, nn, fr;
        if (!num(v[cZ], z) || !num(v[cN], nn) || !num(v[cF], fr)) { ++nbad; continue; }
        const double maf = (fr <= 0.5) ? fr : 1.0 - fr;
        ctx.kz[t][it->second]   = (float) ((fr <= 0.5) ? z : -z);
        ctx.kmac[t][it->second] = (float) (2.0 * nn * maf);
        ++nmatch;
    }
    Rcout << "KVIK summaries " << path << ": " << nread << " variants, " << nmatch
          << " matched to the .bim";
    if (nbad) Rcout << " (" << nbad << " rows skipped: missing or non-numeric values)";
    Rcout << "\n";
    if (nmatch == 0) stop("No variant in " + path + " matched the .bim names -- check the Predictor IDs");
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
// Per-phenotype-variant inputs of the exact binary saddlepoint (spa_used = 5).
// y_i = lo_i + delta b_i: lo is variant specific (it carries the PRS offset), b, pi and
// delta are those of the original phenotype.
struct BinExact {
    bool active;
    GenoMat LoF, PiF;                           // n x P: lo and null P(case)
    std::vector< const std::vector<unsigned char>* > b;   // raw 0/1 per trait
    std::vector<double> delta, ncase;
    std::vector<double> vbar;                   // per trait: delta^2 * mean_i pi_i (1 - pi_i)
    double ec_max;
    double grid_p;                              // lattice convolution when saddlepoint p < grid_p
    BinExact() : active(false), ec_max(10.0), grid_p(0.1) {}
};

struct PhenoStats {
    Eigen::MatrixXd Y; GenoMat Yf;
    std::vector<double> Vp, yty, ycross;
    CoherBin cb;
    BinExact bx;
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
    S.bx = BinExact();
    if (binary && !coher && ctx.bin_exact && (int) ctx.braw.size() == P) {
        S.bx.active = true; S.bx.ec_max = ctx.ec_max;
        S.bx.LoF.resize(n, P); S.bx.PiF.resize(n, P);
        S.bx.b.resize(P); S.bx.delta.resize(P); S.bx.ncase.resize(P); S.bx.vbar.assign(P, 0.0);
        for (int t = 0; t < P; ++t) {
            const double dl = ctx.bdelta[t];
            S.bx.b[t] = &ctx.braw[t]; S.bx.delta[t] = dl; S.bx.ncase[t] = ctx.ncase[t];
            double vb = 0.0;
            for (int i = 0; i < n; ++i) {
                S.bx.LoF(i, t) = (float)(S.Y(i, t) - dl * (double) ctx.braw[t][i]);
                S.bx.PiF(i, t) = ctx.bpi[t][i];
                vb += (double) ctx.bpi[t][i] * (1.0 - (double) ctx.bpi[t][i]);
            }
            S.bx.vbar[t] = dl * dl * vb / n;
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
    std::vector<float> colmaf;         // MAF of every column of V (0 = pseudo-marker);
                                       // lets V be re-weighted to another alpha by a
                                       // pure column rescaling, with no second read
    std::vector<int> colbim;           // .bim index of every CATEGORY column (-1 for
                                       // background / flank columns and pseudo-markers)
    // the ORIGINAL variants of each active category (.bim index and MAF): identical to
    // the columns without collapsing; with collapsing they still include the collapsed
    // ultra-rare variants, so ACAT-V and cat_w see every variant individually
    std::vector< std::vector<int> > cat_obim;
    std::vector< std::vector<float> > cat_omaf;
    // exact binary SPA: minor-allele carriers (analysis-row index) and dosages of every
    // column of V (empty for dense columns), and the column sums of V
    std::vector< std::vector<int> > car;
    std::vector< std::vector<unsigned char> > dos;
    std::vector<double> colsum;
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
    // multi-alpha only: the p-value at each alpha, [cat][trait or pair][alpha].
    // p_spa then holds their ACAT combination.
    std::vector< std::vector< std::vector<double> > > p_alpha;
    // ACAT-V (only with KVIK summaries): per alpha [cat][trait][alpha], the
    // ACAT over alphas of the variance-component p-values (p_vc) and of the
    // ACAT-V p-values (p_acv). p_spa then holds the omnibus of all of them.
    std::vector< std::vector< std::vector<double> > > p_acatv;
    std::vector< std::vector<double> > p_vc, p_acv;
    int n_cond = 0; std::string cond_snps;     // cond_common: variants conditioned on
};

// Read a SNP index range, project covariates, apply alpha; return false on error.
static bool read_block(const ChunkContext& ctx, int lo, int hi, Cell& c, bool collapse = false) {
    c.X = GenoMat(ctx.n_inds, 0);
    if (hi <= lo) return true;
    const bool coll = collapse && ctx.collapse_mac > 0;
    if (!read_cell_idx(ctx.bed, ctx.geno_keep, lo, hi, c,
                       coll ? &ctx.snp_cat : 0, coll ? ctx.collapse_mac : 0, 5,
                       ctx.bin_exact ? ctx.car_max : 0)) return false;
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
    if (!read_block(ctx, a, b, tgt, true)) return false;
    const int m_t = (int) tgt.X.cols();
    if (m_t <= 0) return false;
    const bool coll = !tgt.cat.empty();                 // collapsed: columns != .bim variants

    const int A = ctx.n_cat;
    std::vector< std::vector<int> > cat_cols(A);
    std::vector<int> bg_cols;
    for (int j = 0; j < m_t; ++j) {
        const int cid = coll ? tgt.cat[j] : ctx.snp_cat[a + j];
        if (cid >= 0 && cid < A) cat_cols[cid].push_back(j); else bg_cols.push_back(j);
    }
    // original variants per category (for ACAT-V and the alpha-model weights)
    std::vector< std::vector<int> > ob(A);
    std::vector< std::vector<float> > om(A);
    for (int j = 0; j < b - a; ++j) {
        const int cid = coll ? tgt.ocat[j] : ctx.snp_cat[a + j];
        if (cid < 0 || cid >= A) continue;
        ob[cid].push_back(a + j);
        om[cid].push_back(coll ? tgt.omaf[j] : tgt.maf[j]);
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
    cd.colmaf.assign(cd.K, 0.f);
    cd.colbim.assign(cd.K, -1);
    const bool keep_car = ctx.bin_exact;
    if (keep_car) { cd.car.assign(cd.K, std::vector<int>()); cd.dos.assign(cd.K, std::vector<unsigned char>()); }
    // copy the carrier list of source column j of cell C into column k of V
    auto take_car = [&](const Cell& C, int j, int k) {
        if (!keep_car || j >= (int) C.car.size()) return;
        cd.car[k] = C.car[j]; cd.dos[k] = C.dos[j];
    };
    const double e1 = 1.0 + ctx.alpha;
    int off = 0;
    for (size_t k = 0; k < act.size(); ++k) {
        const int c = act[k]; const std::vector<int>& cc = cat_cols[c];
        for (size_t q = 0; q < cc.size(); ++q) {
            cd.V.col(off + (int) q) = tgt.X.col(cc[q]);
            cd.colmaf[off + (int) q] = tgt.maf[cc[q]];
            cd.colbim[off + (int) q] = coll ? tgt.bim[cc[q]] : a + cc[q];
            take_car(tgt, cc[q], off + (int) q);
        }
        // expected heritability under the alpha model, over the ORIGINAL variants
        double w = 0.0;
        for (size_t q = 0; q < om[c].size(); ++q) {
            const double f = om[c][q];
            w += (f > 0.0 && f < 1.0) ? std::pow(2.0 * f * (1.0 - f), e1) : 0.0;
        }
        cd.cat_obim.push_back(ob[c]); cd.cat_omaf.push_back(om[c]);
        off += (int) cc.size();
        cd.cat_m.push_back((int) cc.size()); cd.cat_id.push_back(c);
        cd.cat_name.push_back(ctx.cat_names[c]); cd.cat_w.push_back(w);
    }
    for (size_t q = 0; q < bg_cols.size(); ++q) {
        take_car(tgt, bg_cols[q], off);
        cd.colmaf[off] = tgt.maf[bg_cols[q]]; cd.V.col(off++) = tgt.X.col(bg_cols[q]);
    }
    if (fl.X.cols() > 0) {
        for (int q = 0; q < (int) fl.X.cols(); ++q) { cd.colmaf[off + q] = fl.maf[q]; take_car(fl, q, off + q); }
        cd.V.middleCols(off, fl.X.cols()) = fl.X; off += (int) fl.X.cols();
    }
    if (fr.X.cols() > 0) {
        for (int q = 0; q < (int) fr.X.cols(); ++q) { cd.colmaf[off + q] = fr.maf[q]; take_car(fr, q, off + q); }
        cd.V.middleCols(off, fr.X.cols()) = fr.X; off += (int) fr.X.cols();
    }
    cd.m_flank = m_flank;
    cd.cj.resize(cd.K);
    for (int j = 0; j < cd.K; ++j) cd.cj[j] = cd.V.col(j).cast<double>().squaredNorm();
    if (keep_car) {
        cd.colsum.resize(cd.K);
        for (int j = 0; j < cd.K; ++j) cd.colsum[j] = cd.V.col(j).cast<double>().sum();
    }
    return true;
}


// ===========================================================================
// CONDITIONING ON THE GENE'S LEAD COMMON VARIANTS (cond_common = TRUE)
// ===========================================================================
// Rare variants sit on haplotypes that carry common alleles, so a gene whose COMMON
// variants are causal shows rare-variant signal through LD, even in a model that also
// contains a common-variant category (a single infinitesimal common kernel cannot
// absorb a few strong causal variants). The standard remedy, as in SAIGE-GENE+ and
// REGENIE, is to condition on the lead common variants: they become fixed effects.
//
// Selection (per gene block): common variants of the block (MAF >= cond_maf) whose
// LDAK-KVIK single-variant p-value is < cond_p for ANY trait, in order of their smallest
// p-value; a variant is skipped if it is (nearly) a linear combination of those already
// selected (r^2 > 0.9 with their span); at most cond_max variants. With Q an orthonormal
// basis of the selected (covariate-projected, weighted) genotype columns:
//   V <- (I - QQ') V      (the selected columns themselves are dropped)
//   y <- (I - QQ') y      (via the offset in build_pheno_stats, every phenotype)
//   T(env, env) = n - cov_df - r
//   binary: P(case) <- P(case) + QQ'(b - P(case))  (linear-probability update), for the
//           exact single-trait test (PiF) and the conditional co-heritability test (cpi).
// Returns r (0 = nothing selected; cd unchanged).
struct CondOpt { double p, maf; int max; };
// p-values: LDAK-KVIK single-variant tests when kvik_summaries are given; otherwise the
// correlation test of each common variant with each phenotype in the data (Y: the
// phenotypes as analysed, i.e. after covariates and any LOCO PRS), z = r sqrt(n).
static int cond_common_apply(const ChunkContext& ctx, const CondOpt& pr, ChunkDataA& cd,
                             Eigen::MatrixXd& Q, std::string& names, const Eigen::MatrixXd* Y = 0) {
    names.clear();
    const bool use_kvik = !ctx.kz.empty();
    const int n = (int) cd.V.rows(), P = use_kvik ? (int) ctx.kz.size() : (Y ? (int) Y->cols() : 0);
    if (P == 0 || n <= 0) return 0;
    Eigen::RowVectorXd ynorm;
    if (!use_kvik) {
        ynorm.resize(P);
        for (int t = 0; t < P; ++t) { const double m = Y->col(t).mean(); ynorm[t] = std::sqrt((Y->col(t).array() - m).square().sum()); }
    }
    int mcat = 0; for (size_t c = 0; c < cd.cat_m.size(); ++c) mcat += cd.cat_m[c];
    std::vector< std::pair<double,int> > cand;
    for (int k = 0; k < mcat; ++k) {
        const int bi = cd.colbim[k];
        if (bi < 0 || !(cd.colmaf[k] >= pr.maf)) continue;
        double pmin = 1.0;
        if (use_kvik) {
            for (int t = 0; t < P; ++t) {
                const float z = ctx.kz[t][bi];
                if (std::isnan(z)) continue;
                pmin = std::min(pmin, std::erfc(std::abs((double) z) / std::sqrt(2.0)));
            }
        } else {
            Eigen::VectorXd x = cd.V.col(k).cast<double>(); x.array() -= x.mean();
            const double xn = x.norm();
            if (!(xn > 0.0)) continue;
            const Eigen::RowVectorXd xy = x.transpose() * (*Y);       // centring x is enough
            for (int t = 0; t < P; ++t) {
                if (!(ynorm[t] > 0.0)) continue;
                const double z = xy[t] / (xn * ynorm[t]) * std::sqrt((double) n);
                pmin = std::min(pmin, std::erfc(std::abs(z) / std::sqrt(2.0)));
            }
        }
        if (pmin < pr.p) cand.push_back(std::make_pair(pmin, k));
    }
    if (cand.empty()) return 0;
    std::sort(cand.begin(), cand.end());
    std::vector<Eigen::VectorXd> qs; std::vector<int> sel;
    for (size_t c = 0; c < cand.size() && (int) sel.size() < pr.max; ++c) {
        Eigen::VectorXd x = cd.V.col(cand[c].second).cast<double>();
        const double x2 = x.squaredNorm();
        if (!(x2 > 0.0)) continue;
        for (int it = 0; it < 2; ++it)                       // Gram-Schmidt, twice for stability
            for (size_t q = 0; q < qs.size(); ++q) x -= qs[q] * qs[q].dot(x);
        const double r2 = x.squaredNorm();
        if (r2 < 0.1 * x2) continue;                         // r^2 > 0.9 with the selected span
        qs.push_back(x / std::sqrt(r2)); sel.push_back(cand[c].second);
    }
    const int r = (int) sel.size();
    if (r == 0) return 0;
    Q.resize(n, r);
    for (int q = 0; q < r; ++q) Q.col(q) = qs[q];
    for (int q = 0; q < r; ++q) { if (q) names += ","; names += ctx.bim.id[cd.colbim[sel[q]]]; }

    // drop the selected columns, project the rest (blocks of 64 columns)
    std::vector<char> drop(cd.K, 0); for (int q = 0; q < r; ++q) drop[sel[q]] = 1;
    std::vector<int> keep; for (int k = 0; k < cd.K; ++k) if (!drop[k]) keep.push_back(k);
    const int K2 = (int) keep.size();
    GenoMat V2(n, K2);
    for (int k0 = 0; k0 < K2; k0 += 64) {
        const int kb = std::min(64, K2 - k0);
        Eigen::MatrixXd B(n, kb);
        for (int c = 0; c < kb; ++c) B.col(c) = cd.V.col(keep[k0 + c]).cast<double>();
        B.noalias() -= Q * (Q.transpose() * B);
        V2.middleCols(k0, kb) = B.cast<float>();
    }
    // per-category bookkeeping
    std::vector<int> cm2; std::vector<int> cid2; std::vector<std::string> cn2; std::vector<double> cw2;
    std::vector< std::vector<int> > cob2; std::vector< std::vector<float> > com2;
    int off = 0;
    for (size_t c = 0; c < cd.cat_m.size(); ++c) {
        int m2 = 0; for (int k = off; k < off + cd.cat_m[c]; ++k) if (!drop[k]) ++m2;
        off += cd.cat_m[c];
        if (m2 == 0) continue;                               // category emptied: no longer tested
        cm2.push_back(m2); cid2.push_back(cd.cat_id[c]); cn2.push_back(cd.cat_name[c]); cw2.push_back(cd.cat_w[c]);
        if (c < cd.cat_obim.size()) { cob2.push_back(cd.cat_obim[c]); com2.push_back(cd.cat_omaf[c]); }
    }
    if (cm2.empty()) return -1;
    std::vector<float> mf2(K2); std::vector<int> cb2(K2);
    std::vector< std::vector<int> > car2; std::vector< std::vector<unsigned char> > dos2;
    const bool hc = !cd.car.empty();
    if (hc) { car2.resize(K2); dos2.resize(K2); }
    for (int k = 0; k < K2; ++k) {
        mf2[k] = cd.colmaf[keep[k]]; cb2[k] = cd.colbim[keep[k]];
        if (hc) { car2[k] = cd.car[keep[k]]; dos2[k] = cd.dos[keep[k]]; }
    }
    cd.V = V2; cd.K = K2; cd.colmaf = mf2; cd.colbim = cb2;
    if (hc) { cd.car = car2; cd.dos = dos2; }
    cd.cat_m = cm2; cd.cat_id = cid2; cd.cat_name = cn2; cd.cat_w = cw2; cd.cat_obim = cob2; cd.cat_omaf = com2;
    cd.cj.resize(K2);
    for (int j = 0; j < K2; ++j) cd.cj[j] = cd.V.col(j).cast<double>().squaredNorm();
    if (!cd.colsum.empty() || hc) {
        cd.colsum.resize(K2);
        for (int j = 0; j < K2; ++j) cd.colsum[j] = cd.V.col(j).cast<double>().sum();
    }
    return r;
}

// Phenotype statistics conditional on the selected variants (basis Q): y <- (I - QQ')y,
// and the binary null probabilities updated by the linear-probability fit on Q.
static void cond_pheno_stats(const ChunkContext& ctx, bool binary, bool coher,
                             const std::vector< std::pair<int,int> >& pairs, const PhenoStats& S0,
                             const Eigen::MatrixXd& Q, PhenoStats& S) {
    const int n = ctx.n_inds, P = ctx.n_pheno;
    Eigen::MatrixXd off = ctx.Y - S0.Y;                       // the LOCO offset already applied
    off.noalias() += Q * (Q.transpose() * S0.Y);
    build_pheno_stats(ctx, off, binary, coher, pairs, S);
    auto upd = [&](std::vector<double>& pv, const std::vector<unsigned char>& b) {
        Eigen::VectorXd d(n);
        for (int i = 0; i < n; ++i) d[i] = (double) b[i] - pv[i];
        const Eigen::VectorXd a = Q * (Q.transpose() * d);
        for (int i = 0; i < n; ++i) pv[i] = std::min(std::max(pv[i] + a[i], 1e-8), 1.0 - 1e-8);
    };
    if (S.bx.active) {
        for (int t = 0; t < P; ++t) {
            std::vector<double> pv(n);
            for (int i = 0; i < n; ++i) pv[i] = S.bx.PiF(i, t);
            upd(pv, ctx.braw[t]);
            double vb = 0.0;
            for (int i = 0; i < n; ++i) { S.bx.PiF(i, t) = (float) pv[i]; vb += pv[i] * (1.0 - pv[i]); }
            S.bx.vbar[t] = S.bx.delta[t] * S.bx.delta[t] * vb / n;
        }
    }
    if (S.cb.active) {
        for (size_t q = 0; q < S.cb.cpi.size(); ++q) {
            std::vector<double> pv(n);
            for (int i = 0; i < n; ++i) pv[i] = S.cb.cpi[q][i];
            upd(pv, ctx.braw[S.cb.rand_t[q]]);
            for (int i = 0; i < n; ++i) S.cb.cpi[q][i] = (float) pv[i];
        }
    }
}

// ===========================================================================
// Chunk geometry: everything that depends on the GENOTYPES and cov_df but NOT
// on the phenotype -- the Gram, G^2, the moment matrix T and its decomposition,
// and (when SPA is on) the Cholesky factor L and the block outer products W.
//
// This is ~97% of the cost of a window at K ~ 3000: the Gram alone is nK^2/2
// flops against nKP for q = V'y. Building it once and sweeping the phenotype
// variants (distinct prs_mask rows) over it is what makes pass 2 cost ONE
// traversal of the window list instead of one per variant. Results are
// unchanged: the same operations in the same order, just not repeated.
// ===========================================================================
struct ChunkGeom {
    bool ok;
    int K, A, C, env;
    bool has_f;
    std::vector<int> off;            // component column offsets in V
    std::vector<double> g;           // trace normalisers n / tr(S_a)
    Eigen::MatrixXd Tinv;
    Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> Tcod;
    // SPA scaffolding -- also phenotype-independent
    bool have_L, use_wcache;
    Eigen::MatrixXd L;
    std::vector<Eigen::MatrixXd> W;
    ChunkGeom() : ok(false), K(0), A(0), C(0), env(0), has_f(false),
                  have_L(false), use_wcache(false) {}
};

// Gpre: an already-formed (double) Gram for this cd, used instead of recomputing
// V'V. Gout: if non-null, receives the Gram. Both default to null, in which case
// this is exactly the original function. They exist so the multi-alpha path can
// form the Gram ONCE and obtain every other alpha's as  D G D  -- a K^2 rescale
// instead of another n K^2 / 2 product.
static void build_geom(const ChunkDataA& cd, int cov_df, bool spa, bool binary,
                       ChunkGeom& gm, const Eigen::MatrixXd* Gpre = 0,
                       Eigen::MatrixXd* Gout = 0) {
    gm = ChunkGeom();
    const int n = (int) cd.V.rows(), K = cd.K;
    const int A = (int) cd.cat_m.size();
    if (A <= 0 || K <= 0 || n <= 0) return;
    const bool has_f = (cd.m_flank > 0);
    // Components: the A tested categories, then the flank ONLY if present. With no
    // flank (pass 1) the model is { cat_0, ..., cat_{A-1}, sigma_e I }.
    const int C = A + (has_f ? 1 : 0), env = C;
    gm.K = K; gm.A = A; gm.C = C; gm.env = env; gm.has_f = has_f;

    gm.off.assign(C + 1, 0);
    for (int c = 0; c < A; ++c) gm.off[c + 1] = gm.off[c] + cd.cat_m[c];
    if (has_f) gm.off[A + 1] = gm.off[A] + cd.m_flank;

    Eigen::MatrixXd G;
    if (Gpre != 0 && Gpre->rows() == K && Gpre->cols() == K) {
        G = *Gpre;
    } else {
        GenoMat Gf = GenoMat::Zero(K, K);
        Gf.selfadjointView<Eigen::Upper>().rankUpdate(cd.V.transpose());
        G = GenoMat(Gf.selfadjointView<Eigen::Upper>()).cast<double>();
    }
    if (Gout != 0) *Gout = G;

    {
        Eigen::MatrixXd G2 = G.array().square();
        std::vector<double> tr(C, 0.0);
        gm.g.assign(C, 0.0);
        for (int c = 0; c < C; ++c) {
            for (int j = gm.off[c]; j < gm.off[c + 1]; ++j) tr[c] += cd.cj[j];
            if (!(tr[c] > 0.0)) return;          // gm.ok stays false
            gm.g[c] = (double) n / tr[c];
        }
        Eigen::MatrixXd T = Eigen::MatrixXd::Zero(C + 1, C + 1);
        for (int a = 0; a < C; ++a)
            for (int b = a; b < C; ++b) {
                double v = gm.g[a] * gm.g[b] *
                    G2.block(gm.off[a], gm.off[b],
                             gm.off[a+1] - gm.off[a], gm.off[b+1] - gm.off[b]).sum();
                T(a, b) = v; T(b, a) = v;
            }
        // See test_chunk: tr(P) = n - cov_df, not n, once covariates are regressed out.
        for (int a = 0; a < C; ++a) { T(a, env) = (double) n; T(env, a) = (double) n; }
        T(env, env) = (double) (n - cov_df);
        gm.Tcod.compute(T);
        gm.Tinv = gm.Tcod.pseudoInverse();
    }

    if (spa) {
        double gscale = G.trace() / K;
        Eigen::MatrixXd Gj = G; double applied = 0.0;
        for (double r = 1e-4; r <= 1e-1; r *= 10.0) {
            Gj.diagonal().array() += (r - applied) * gscale; applied = r;
            Eigen::LLT<Eigen::MatrixXd> llt(Gj);
            if (llt.info() == Eigen::Success) { gm.L = llt.matrixL(); gm.have_L = true; break; }
        }
    }

    // ---- block outer products, once per chunk (shared across traits AND
    // categories -- the saving is A x P fold here, not just P fold) ----------
    gm.use_wcache = spa && gm.have_L && ((double) C * K * K * 8.0 <= 6e8);
    if (gm.use_wcache) {
        gm.W.resize(C);
        for (int a = 0; a < C; ++a) {
            const int ma = gm.off[a + 1] - gm.off[a];
            gm.W[a].noalias() = gm.L.middleRows(gm.off[a], ma).transpose()
                              * gm.L.middleRows(gm.off[a], ma);
        }
    }

    gm.ok = true;
}


static void test_chunk_annot(const ChunkDataA& cd, ChunkGeom& gm,
                             const Eigen::MatrixXd& Y, const GenoMat& Yf,
                             const std::vector<double>& Vp, const std::vector<double>& yty,
                             bool binary,
                             bool spa, bool coher,
                             const std::vector< std::pair<int,int> >& pairs,
                             const std::vector<double>& ycross,
                             const CoherBin& cb,
                             const BinExact& bx, bool one_sided,
                             ChunkResultA& cr) {
    const int n = (int) Y.rows(), P = (int) Y.cols();
    const int A = (int) cd.cat_m.size();

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
    if (!gm.ok) return;          // degenerate geometry: metadata + NAs only

    // REPORTED STANDARD ERROR. The closed-form variance of the moment estimator
    // (var0 below) is still computed -- the saddlepoint needs it internally to
    // normalise the null distribution and pin its second cumulant -- but it is
    // NOT reported: for rare-variant kernels it does not describe the sampling
    // distribution of vg well. The reported se_vg is instead the SE implied by
    // the estimate and its saddlepoint p-value, se = |vg| / z(p), so that
    // vg / se_vg reproduces the p-value exactly.
    auto finish_se = [&]() {
        for (size_t c = 0; c < cr.vg.size(); ++c)
            for (size_t t = 0; t < cr.vg[c].size(); ++t)
            {
                // the SE always comes from the TWO-SIDED p (defined for negative estimates too)
                const double pp = cr.p_spa[c][t];
                const double p2 = (one_sided && !coher && std::isfinite(pp)) ? std::min(1.0, 2.0 * std::min(pp, 1.0 - pp)) : pp;
                cr.se_vg[c][t] = se_from_p(cr.vg[c][t], p2);
            }
    };

    // Aliases onto the prebuilt geometry. Everything below is exactly as it was
    // when these were local; they are simply not recomputed per variant.
    const int K = gm.K, C = gm.C, env = gm.env;
    const bool has_f = gm.has_f, have_L = gm.have_L, use_wcache = gm.use_wcache;
    const std::vector<int>& off = gm.off;
    const std::vector<double>& g = gm.g;
    const Eigen::MatrixXd& Tinv = gm.Tinv;
    const Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd>& Tcod = gm.Tcod;
    const Eigen::MatrixXd& L = gm.L;
    const std::vector<Eigen::MatrixXd>& W = gm.W;

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


    Eigen::MatrixXd U = (cd.V.transpose() * Yf).cast<double>();      // K x P, one GEMM

    // ---- exact binary SPA (spa_used = 5): shifts and Gaussian means, one GEMM each ----
    const bool use_exact = binary && !coher && spa && have_L && bx.active &&
                           (int) cd.car.size() == K && (int) cd.colsum.size() == K;
    Eigen::MatrixXd Ulo, Upi;
    std::vector<int> comp_of(K, env);
    if (use_exact) {
        Ulo = (cd.V.transpose() * bx.LoF).cast<double>();
        Upi = (cd.V.transpose() * bx.PiF).cast<double>();
        for (int a = 0; a < C; ++a) for (int j = off[a]; j < off[a + 1]; ++j) comp_of[j] = a;
    }

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
        // ---- conditional null (see the supplement): condition on the first trait y1 ----
        // sigma12_hat = y1' A_c y2 = h' y2 with h = A_c y1, which given y1 is LINEAR in y2.
        // Per-trait quantities that do not depend on the pair are cached:
        //   a_t = L^{-1} X' y_t (coordinates of y_t in the column space of X), ||r_t||^2 the
        //   squared norm of its orthogonal part, and the Cholesky factor of Psi_t (trait t's
        //   m-dimensional covariance part), and, for binary traits, h = A_c y_t (an n-vector).
        std::vector<Eigen::VectorXd> acache(P); std::vector<double> rcache(P, 0.0), e1cache(P, 1.0);
        std::vector<char> ahave(P, 0), phave(P, 0);
        std::vector< Eigen::LLT<Eigen::MatrixXd> > pcache(P);
        auto prep_trait = [&](int t) {
            if (ahave[t]) return;
            acache[t] = L.triangularView<Eigen::Lower>().solve(U.col(t));
            rcache[t] = std::max(0.0, yty[t] - acache[t].squaredNorm());
            ahave[t] = 1;
        };
        // Psi1 uses the UNCLAMPED estimates of trait 1 (clamping at 0 biases the conditional
        // mean in one direction when the true component is 0; the estimation error itself is
        // propagated in the test, see coher_reml_test). Clamped only if Psi1 is not PD.
        std::vector<Eigen::MatrixXd> psi1m(P);
        auto prep_psi1 = [&](int t) {
            if (phave[t]) return;
            std::vector<double> d1(C);
            for (int a2 = 0; a2 < C; ++a2) d1[a2] = sg[t][a2] * g[a2];
            e1cache[t] = std::max(sg[t][env], 1e-8 * (Vp[t] > 0 ? Vp[t] : 1.0));
            coher_block(L, W, use_wcache, off, C, K, d1, e1cache[t], psi1m[t]);
            pcache[t].compute(psi1m[t]);
            if (pcache[t].info() != Eigen::Success) {
                for (int a2 = 0; a2 < C; ++a2) d1[a2] = std::max(sg[t][a2], 0.0) * g[a2];
                coher_block(L, W, use_wcache, off, C, K, d1, e1cache[t], psi1m[t]);
                pcache[t].compute(psi1m[t]);
            }
            phave[t] = 1;
        };
        // sampling covariance of trait 1's component estimates (categories 0..A-1):
        //   Cov(sigma_a, sigma_b) = 2 tr(A_a Psi1 A_b Psi1) + 2 (n - K) T(a,env) T(b,env) e1^2
        std::vector<Eigen::MatrixXd> vth(P); std::vector<char> vhave(P, 0);
        auto prep_vth = [&](int t) {
            if (vhave[t]) return;
            prep_psi1(t);
            std::vector<Eigen::MatrixXd> AP(A);
            for (int a2 = 0; a2 < A; ++a2) AP[a2] = Acat[a2] * psi1m[t];
            vth[t].setZero(A, A);
            for (int a2 = 0; a2 < A; ++a2) for (int b2 = a2; b2 < A; ++b2) {
                const double v = 2.0 * (AP[a2].cwiseProduct(AP[b2].transpose())).sum()
                               + 2.0 * (double) (n - K) * Tinv(a2, env) * Tinv(b2, env) * e1cache[t] * e1cache[t];
                vth[t](a2, b2) = vth[t](b2, a2) = v;
            }
            vhave[t] = 1;
        };
        std::vector< std::map<int, Eigen::VectorXd> > hcache(A);
        // ---- small genes: REML-GLS conditional test ------------------------------------
        // Given y1, y2 - mu ~ N(0, Omega(s)) on the column space of X (coordinates a2 - m)
        // and N(0, v_perp I) on its complement, with Omega(s) = P0 - Psi12 Psi1^{-1} Psi12 + s M_c,
        // s the (unknown) tested component of trait 2. A genetic covariance in component c
        // adds a term proportional to h to the conditional mean, so we test the fixed effect
        // of the covariate h (column-space part b, orthogonal part eps_c r1) by generalised
        // least squares, with s estimated by REML (which accounts for the fixed effect, so the
        // estimate is not driven by the h-direction of y2) and a Satterthwaite t for its
        // uncertainty. Omega is diagonalised once per pair (generalised eigenproblem).
        const int KREML = 512;            // = default max_gene_snps: every gene block uses this test
        auto coher_reml_test = [&](int c, int t1, int t2, const Eigen::VectorXd& b, const Eigen::VectorXd& a1,
                                   double r1, double e1, double e2, double e12, const Eigen::VectorXd& w,
                                   const Eigen::MatrixXd& Psi12, const Eigen::MatrixXd& P0, const Eigen::MatrixXd& M,
                                   const Eigen::LLT<Eigen::MatrixXd>& l1, double ec,
                                   const std::vector<double>& yty_, double y12,
                                   const std::vector<Eigen::VectorXd>& msl, const Eigen::MatrixXd& Vth,
                                   const std::vector<Eigen::VectorXd>& xbg, double& tst, double& nu) -> bool {
            (void) t1; (void) w; (void) yty_;
            prep_trait(t2);
            const Eigen::VectorXd& a2 = acache[t2];
            const double vperp = e2 - e12 * e12 / e1;
            if (!(vperp > 0.0)) return false;
            // conditional mean and covariance on the column space
            const Eigen::VectorXd z = a2 - Psi12 * l1.solve(a1);
            Eigen::MatrixXd P0c = P0 - Psi12 * l1.solve(Psi12);
            P0c = (0.5 * (P0c + P0c.transpose())).eval();
            Eigen::LLT<Eigen::MatrixXd> lp(P0c);
            if (lp.info() != Eigen::Success) return false;
            const Eigen::MatrixXd Fi = lp.matrixL().solve(Eigen::MatrixXd::Identity(K, K));   // F^{-1}
            Eigen::MatrixXd Sm = Fi * M * Fi.transpose(); Sm = (0.5 * (Sm + Sm.transpose())).eval();
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Sm);
            if (es.info() != Eigen::Success) return false;
            const Eigen::VectorXd lam = es.eigenvalues().cwiseMax(0.0);
            const Eigen::MatrixXd Vt = es.eigenvectors().transpose() * Fi;      // whitening + eigenbasis
            const Eigen::VectorXd zt = Vt * z;
            // FIXED EFFECTS. Column 0 = the tested covariate h. The cross-covariances that the
            // conditional mean was built from (e12 and the genetic covariances of the other
            // components of the block) are estimated from the same two phenotypes; plugging them
            // in ignores their estimation error, which inflates the test when the phenotypes are
            // correlated. They are therefore re-estimated JOINTLY with the test, as nuisance
            // covariates in the same GLS (z is the residual after the plug-in mean, so their
            // coefficients estimate the plug-in errors):
            //   env:          column part Psi1^{-1} a1, orthogonal part |r1|^2 / e1
            //   component a:  column part B_a Psi1^{-1} a1  (B_a = kernel of component a), a != c
            // Covariates are only added while at least 4 dimensions remain for the residual.
            const Eigen::VectorXd p1a = l1.solve(a1);
            std::vector<Eigen::VectorXd> xc; std::vector<double> xo;    // column / orthogonal parts
            const double sr1 = std::sqrt(std::max(r1, 0.0));
            xc.push_back(b); xo.push_back(ec * sr1);
            xc.push_back(p1a); xo.push_back(sr1 / e1);
            for (size_t q = 0; q < xbg.size(); ++q) {
                if ((int) xc.size() + 4 > K + 1) break;
                xc.push_back(xbg[q]); xo.push_back(0.0);
            }
            const int NX = (int) xc.size();
            Eigen::MatrixXd Xt(K, NX); Eigen::VectorXd Xo(NX);
            for (int q = 0; q < NX; ++q) { Xt.col(q) = Vt * xc[q]; Xo[q] = xo[q]; }
            // UNCERTAINTY OF TRAIT 1's COMPONENTS. The conditional mean Psi12 Psi1^{-1} a1 uses the
            // components of trait 1 estimated from y1; their errors move the mean (for the tested
            // component almost exactly along the tested direction, when the phenotypes are
            // correlated), which no covariate can absorb without absorbing the signal. Their effect
            // is propagated by the delta method:
            //   beta_hat changes by -sum_a (d beta / d z)' m_a (sigma1_a_hat - sigma1_a),
            //   m_a = d mean / d sigma1_a = -Psi12 Psi1^{-1} (g_a B_a) Psi1^{-1} a1,
            // which adds  s' Cov(sigma1_hat) s  (s_a = (d beta / d z)' m_a) to Var(beta_hat).
            const int NA = (int) msl.size();
            std::vector<Eigen::VectorXd> mtl(NA);
            for (int q = 0; q < NA; ++q) mtl[q] = Vt * msl[q];
            // orthogonal part: response r1'(y2 - mu) / |r1|, variance v_perp
            const double zo = ((y12 - a1.dot(a2)) - (e12 / e1) * r1) / std::max(sr1, 1e-300);
            auto parts = [&](double sv, Eigen::VectorXd& N, Eigen::MatrixXd& D, double& zz, double& ld) {
                N = Xo * (zo / vperp); D = Xo * Xo.transpose() / vperp; zz = zo * zo / vperp; ld = 0.0;
                for (int k = 0; k < K; ++k) {
                    const double d = 1.0 + sv * lam[k];
                    N.noalias() += Xt.row(k).transpose() * (zt[k] / d);
                    D.noalias() += Xt.row(k).transpose() * Xt.row(k) / d;
                    zz += zt[k] * zt[k] / d; ld += std::log(d);
                }
            };
            auto gls = [&](double sv, double& beta, double& var, double& rl) -> bool {   // beta_h, Var, REML
                Eigen::VectorXd N; Eigen::MatrixXd D; double zz, ld; parts(sv, N, D, zz, ld);
                Eigen::LDLT<Eigen::MatrixXd> ldl(D);
                if (ldl.info() != Eigen::Success || !(ldl.vectorD().minCoeff() > 0.0)) return false;
                const Eigen::VectorXd bh = ldl.solve(N);
                Eigen::VectorXd e0 = Eigen::VectorXd::Zero(NX); e0[0] = 1.0;
                const Eigen::VectorXd de = ldl.solve(e0);
                var = de[0]; beta = bh[0];
                if (NA > 0) {
                    Eigen::VectorXd sens = Eigen::VectorXd::Zero(NA);
                    for (int k = 0; k < K; ++k) {
                        const double wk = Xt.row(k).dot(de) / (1.0 + sv * lam[k]);
                        for (int q = 0; q < NA; ++q) sens[q] += wk * mtl[q][k];
                    }
                    var += sens.dot(Vth * sens);
                }
                rl = -0.5 * (ld + zz - N.dot(bh) + ldl.vectorD().array().log().sum());
                return var > 0.0 && std::isfinite(beta);
            };
            auto reml = [&](double sv) { double bb, vv, rl; return gls(sv, bb, vv, rl) ? rl : -1e300; };
            auto varb = [&](double sv) { double bb, vv, rl; return gls(sv, bb, vv, rl) ? vv : NA_REAL; };
            // maximise over s >= 0 on a log grid, then golden-section refinement
            const double lmax = lam.maxCoeff();
            if (!(lmax > 0.0)) {                              // no information on s: plain GLS
                double bb, vv, rl; if (!gls(0.0, bb, vv, rl)) return false;
                tst = bb / std::sqrt(vv); nu = 1e300; return std::isfinite(tst);
            }
            double best = 0.0, fbest = reml(0.0);
            for (int k = -6; k <= 6; ++k) {
                const double sv = std::pow(10.0, k) / lmax, f = reml(sv);
                if (f > fbest) { fbest = f; best = sv; }
            }
            double lo = (best > 0.0) ? best / 10.0 : 0.0, hi = (best > 0.0) ? best * 10.0 : 1e-6 / lmax;
            const double gr = 0.5 * (std::sqrt(5.0) - 1.0);
            double x1 = hi - gr * (hi - lo), x2 = lo + gr * (hi - lo), f1 = reml(x1), f2 = reml(x2);
            for (int it = 0; it < 60; ++it) {
                if (f1 > f2) { hi = x2; x2 = x1; f2 = f1; x1 = hi - gr * (hi - lo); f1 = reml(x1); }
                else         { lo = x1; x1 = x2; f1 = f2; x2 = lo + gr * (hi - lo); f2 = reml(x2); }
            }
            double sh = (f1 > f2) ? x1 : x2;
            if (reml(0.0) >= std::max(f1, f2)) sh = 0.0;
            double bh, vh, rlh;
            if (!gls(sh, bh, vh, rlh)) return false;
            tst = bh / std::sqrt(vh);
            // Satterthwaite degrees of freedom for the uncertainty of s:
            //   nu = 2 Var(beta_hat)^2 info(s) / (d Var(beta_hat) / ds)^2   at s = s_hat.
            // At (or near) the boundary s_hat = 0 this breaks down. In SMALL blocks (K <= 30), s_hat = 0
            // can badly underestimate a real heritability of trait 2 at the gene (inflation up to ~8x
            // at p = 1e-3 in simulations with h2 = 0.0002 for a 6-variant gene), so there we use
            // nu = K - 1 (conservative for small genes without heritability, which costs little power:
            // co-heritability needs both traits to be heritable at the gene). In larger blocks s is
            // estimated precisely enough and the test at s_hat = 0 was calibrated (and was very
            // conservative with nu = K - 1), so it is normal.
            auto curv = [&](double sv) {                       // -d2 REML / ds2 at sv (one-sided at 0)
                const double hs = std::max(1e-3 * sv, 1e-4 / lmax);
                const double f0 = reml(sv), fp = reml(sv + hs);
                if (sv - hs >= 0.0) return -(fp - 2.0 * f0 + reml(sv - hs)) / (hs * hs);
                return -(reml(sv + 2.0 * hs) - 2.0 * fp + f0) / (hs * hs);
            };
            const double info0 = curv(0.0);
            const double se0 = (info0 > 0.0) ? 1.0 / std::sqrt(info0) : 1.0 / lmax;
            if (sh > 0.1 * se0) {
                nu = std::max(1.0, (double) K - 1.0);
                const double info = curv(sh);
                const double hs = std::max(1e-3 * sh, 1e-4 / lmax);
                const double v1 = varb(sh), v2 = varb(sh + hs), dv = (v2 - v1) / hs;   // d Var(beta_hat) / ds
                if (info > 0.0 && std::isfinite(dv) && std::abs(dv) > 0.0) nu = 2.0 * v1 * v1 * info / (dv * dv);
            } else {
                nu = (K <= 30) ? std::max(1.0, (double) K - 1.0) : 1e300;
            }
            if (nu < 1.0) nu = 1.0;
            return std::isfinite(tst);
        };
        // trait-2 parts per (trait, category): P0 = m-dimensional covariance WITHOUT the tested
        // component (background components and e2 from y2 itself), and the traces giving
        // Var(sigma2_c hat) = 2 (tau0 + tau1 s + tau2 s^2) as a function of the tested
        // component s (univariate null variance 2 tr((A_c Sigma2)^2), expanded in s)
        std::vector<Eigen::MatrixXd> Mc(A);
        std::vector<double> tau2c(A, 0.0);
        std::vector<char> mhave(A, 0);
        struct T2 { Eigen::MatrixXd P0; double e2, tau0, tau1; };
        std::vector< std::map<int, T2> > t2cache(A);
        auto get_M = [&](int c) {
            if (mhave[c]) return;
            std::vector<double> dt(C, 0.0); dt[c] = g[c];
            coher_block(L, W, use_wcache, off, C, K, dt, 0.0, Mc[c]);
            Eigen::MatrixXd RM = Acat[c] * Mc[c];
            tau2c[c] = (RM.cwiseProduct(RM.transpose())).sum();
            mhave[c] = 1;
        };
        auto get_t2 = [&](int c, int t) -> const T2& {
            std::map<int, T2>::iterator it = t2cache[c].find(t);
            if (it == t2cache[c].end()) {
                get_M(c);
                T2 v;
                std::vector<double> d2(C);
                for (int a2 = 0; a2 < C; ++a2) d2[a2] = (a2 == c) ? 0.0 : std::max(sg[t][a2], 0.0) * g[a2];
                v.e2 = std::max(sg[t][env], 1e-8 * (Vp[t] > 0 ? Vp[t] : 1.0));
                coher_block(L, W, use_wcache, off, C, K, d2, v.e2, v.P0);
                Eigen::MatrixXd RP = Acat[c] * v.P0, RM = Acat[c] * Mc[c];
                const double ec = Tinv(c, env);
                v.tau0 = (RP.cwiseProduct(RP.transpose())).sum() + (double) (n - K) * ec * ec * v.e2 * v.e2;
                v.tau1 = 2.0 * (RP.cwiseProduct(RM.transpose())).sum();
                it = t2cache[c].insert(std::make_pair(t, v)).first;
            }
            return it->second;
        };
        auto get_h = [&](int c, int t) -> const Eigen::VectorXd& {
            std::map<int, Eigen::VectorXd>::iterator it = hcache[c].find(t);
            if (it == hcache[c].end()) {
                Eigen::VectorXd wv;
                coher_wvec(cd.V, off, C, K, g, Tinv, c, env, U.col(t), Y.col(t), wv);
                it = hcache[c].insert(std::make_pair(t, wv)).first;
            }
            return it->second;
        };
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
                const double Q = s12[c];
                // ---- BINARY: exact conditional distribution given the lower-prevalence trait ----
                // Q = c0 + sum_i delta h_i b_i, b_i ~ Bernoulli(pi_i) independently. The tail is
                // driven by the few individuals with large |h_i| (carriers): up to 500 of them,
                // with |delta h_i| >= 0.05 SD, are kept as exact two-point terms; the rest enters
                // through its exact (grouped) Bernoulli CGF. The terms are combined by saddlepoint
                // or, when that fails or p < 0.1, by the lattice convolution of coher_grid_solve.
                if (binary && cb.active) {
                    const int tc = cb.cond_t[pi];
                    const Eigen::VectorXd& wv = get_h(c, tc);
                    const double dl = cb.cdelta[pi];
                    const std::vector<float>& cy = cb.cylo[pi];
                    const std::vector<float>& cp = cb.cpi[pi];
                    double c0 = 0.0, sd2 = 0.0;
                    for (int i = 0; i < n; ++i) {
                        c0 += wv[i] * cy[i];
                        const double ci = dl * wv[i];
                        sd2 += ci * ci * cp[i] * (1.0 - cp[i]);
                    }
                    const double thr = 0.05 * std::sqrt(sd2);
                    std::vector<int> big;
                    for (int i = 0; i < n; ++i) if (std::abs(dl * wv[i]) > thr) big.push_back(i);
                    const size_t MAXT = 500;
                    if (big.size() > MAXT) {
                        std::nth_element(big.begin(), big.begin() + MAXT, big.end(),
                                         [&](int x, int y) { return std::abs(wv[x]) > std::abs(wv[y]); });
                        big.resize(MAXT);
                    }
                    std::vector<char> isbig(n, 0);
                    std::vector<ExactTerm> terms; terms.reserve(big.size());
                    for (size_t r = 0; r < big.size(); ++r) {
                        const int i = big[r]; isbig[i] = 1;
                        const double pr = std::min(std::max((double) cp[i], 1e-12), 1.0 - 1e-12);
                        ExactTerm e; e.v.resize(2); e.lp.resize(2);
                        e.v[0] = 0.0; e.lp[0] = std::log1p(-pr);
                        e.v[1] = dl * wv[i]; e.lp[1] = std::log(pr);
                        terms.push_back(e);
                    }
                    // the rest: EXACT CGF of its Bernoulli terms, grouped by (pi, coefficient rounded
                    // to 1e-3 of its SD; the group uses the mean coefficient). A Gaussian for the rest
                    // is not good enough at low prevalence: the rest is then a sum of many
                    // Bernoulli(0.01) terms, i.e. Poisson-like and skewed (anti-conservative tails).
                    double vr = 0.0;
                    for (int i = 0; i < n; ++i) {
                        if (isbig[i]) continue;
                        const double ci = dl * wv[i];
                        vr += ci * ci * cp[i] * (1.0 - cp[i]);
                    }
                    BernRest br;
                    if (vr > 0.0) {
                        const double step = 1e-3 * std::sqrt(vr);
                        std::map< std::pair<float, long>, std::pair<double, double> > grp;
                        for (int i = 0; i < n; ++i) {
                            if (isbig[i]) continue;
                            const double ci = dl * wv[i];
                            if (ci == 0.0) continue;
                            std::pair<double, double>& gg = grp[std::make_pair(cp[i], std::lround(ci / step))];
                            gg.first += ci; gg.second += 1.0;
                        }
                        br.c.reserve(grp.size()); br.p.reserve(grp.size()); br.w.reserve(grp.size());
                        for (std::map< std::pair<float, long>, std::pair<double, double> >::const_iterator it = grp.begin(); it != grp.end(); ++it) {
                            br.c.push_back(it->second.first / it->second.second);
                            br.p.push_back(std::min(std::max((double) it->first.first, 1e-12), 1.0 - 1e-12));
                            br.w.push_back(it->second.second);
                        }
                    }
                    const std::vector<double> nolam;
                    ExactSpaResult er; er.converged = false;
                    if (!terms.empty() || !br.c.empty()) {
                        er = exact_spa_solve(Q, terms, nolam, nolam, c0, &br);
                        // saddlepoint first; for p < 0.1 (or if it fails) the lattice convolution
                        // (coher_grid_solve: terms centred at their most likely value)
                        if (!terms.empty() && !er.exact &&
                            (!er.converged || std::min(er.p_up, er.p_lo) < 0.1)) {
                            std::vector<double> tv(terms.size()), tp(terms.size());
                            for (size_t r = 0; r < terms.size(); ++r) { tv[r] = terms[r].v[1]; tp[r] = std::exp(terms[r].lp[1]); }
                            ExactSpaResult eg = coher_grid_solve(Q, tv, tp, &br, c0);
                            if (eg.converged) er = eg;
                        }
                    }
                    // diagnostics: SPARROW_COHER_DUMP=<dir> writes the null model of every binary
                    // co-heritability test with p < 1e-6 (first 200) for checking against Monte Carlo
                    if (er.converged && 2.0 * std::min(er.p_up, er.p_lo) < 1e-6) {
                        static const char* ddir = std::getenv("SPARROW_COHER_DUMP");
                        static std::atomic<int> ndump(0);
                        if (ddir) { const int id = ndump++; if (id < 200) {
                            char fn[1024]; std::snprintf(fn, sizeof(fn), "%s/coher_dump_%d.txt", ddir, id);
                            if (FILE* f = std::fopen(fn, "w")) {
                                std::fprintf(f, "%.17g %.17g %zu %zu %d %d %d %d %.6g\n", Q, c0, terms.size(), br.c.size(), t1, t2, c, K, 2.0 * std::min(er.p_up, er.p_lo));
                                for (size_t r = 0; r < terms.size(); ++r) std::fprintf(f, "%.17g %.17g\n", terms[r].v[1], std::exp(terms[r].lp[1]));
                                for (size_t g2 = 0; g2 < br.c.size(); ++g2) std::fprintf(f, "%.17g %.17g %.17g\n", br.c[g2], br.p[g2], br.w[g2]);
                                std::fclose(f);
                            } } }
                    }
                    if (er.converged) {
                        cr.p_spa[c][pi] = std::min(1.0, 2.0 * std::min(er.p_up, er.p_lo));
                        cr.spa_used[c][pi] = 4;
                    } else { cr.p_spa[c][pi] = NA_REAL; cr.spa_used[c][pi] = -1; }
                    continue;
                }

                // ---- QUANTITATIVE: conditional null given y1 ----------------------------------
                // y2 | y1 ~ N(S21 S1^{-1} y1, S2 - S21 S1^{-1} S12), so Q | y1 is EXACTLY normal
                // with mean h'mu and variance h'V h. With b = R_c a1 (R_c the m-dimensional part of
                // A_c) and r1 the part of y1 orthogonal to the columns of X (Lemma 1):
                //   h'mu  = b' Psi12 Psi1^{-1} a1 + eps_c (e12 / e1) |r1|^2
                //   h'V h = b' Psi2 b - (Psi12 b)' Psi1^{-1} (Psi12 b) + eps_c^2 |r1|^2 (e2 - e12^2 / e1)
                // Psi1 is a function of y1 (conditioned on). The only nuisance from y2 that matters is
                // the tested component of trait 2 (s), which enters h'V h as s * b' M_c b. It is
                // estimated from y2 with its h-direction removed (so (nearly) independent of h'y2),
                // and its uncertainty is propagated by a Satterthwaite approximation:
                //   (Q - h'mu) / sqrt(V_hat) ~ t_nu,   nu = 2 V_hat^2 / ((b' M_c b)^2 Var(s_hat)).
                // For genes with many variants nu is large and this is the normal test; for small
                // genes where trait 2 is heritable, s_hat has few degrees of freedom and the t
                // distribution prevents the inflation a plug-in normal test would have.
                prep_trait(t1); prep_psi1(t1);
                if (pcache[t1].info() != Eigen::Success) { cr.spa_used[c][pi] = -1; continue; }
                const Eigen::VectorXd& a1 = acache[t1];
                const double r1 = rcache[t1], e1 = e1cache[t1], ec = Tinv(c, env);
                const Eigen::VectorXd b = Acat[c] * a1;
                const double hh = b.squaredNorm() + ec * ec * r1;
                if (!(hh > 0.0)) { cr.spa_used[c][pi] = -1; continue; }
                const T2& t2p = get_t2(c, t2);
                // tested component of trait 2 from y2 with the h-direction removed
                const Eigen::VectorXd u2p = U.col(t2) - (L * b) * (Q / hh);
                Eigen::VectorXd qp(C + 1);
                for (int cc = 0; cc < C; ++cc) {
                    double qa = 0.0;
                    for (int j = off[cc]; j < off[cc + 1]; ++j) qa += u2p[j] * u2p[j];
                    qp[cc] = g[cc] * qa;
                }
                qp[env] = yty[t2] - Q * Q / hh;
                const double s2 = Tcod.solve(qp)[c];
                const double e2 = t2p.e2;
                double e12 = s12[env];
                const double ecap = 0.999 * std::sqrt(e1 * e2);
                if (e12 > ecap) e12 = ecap; if (e12 < -ecap) e12 = -ecap;
                std::vector<double> d12(C);
                for (int a2 = 0; a2 < C; ++a2) d12[a2] = (a2 == c) ? 0.0 : s12[a2] * g[a2];   // null: no covariance in c
                Eigen::MatrixXd Psi12;
                coher_block(L, W, use_wcache, off, C, K, d12, e12, Psi12);
                const Eigen::VectorXd w = Psi12 * b;
                const double hm = w.dot(pcache[t1].solve(a1)) + ec * (e12 / e1) * r1;
                const double V0 = b.dot(t2p.P0 * b) - w.dot(pcache[t1].solve(w)) + ec * ec * r1 * (e2 - e12 * e12 / e1);
                const double bm = b.dot(Mc[c] * b);
                if (!(V0 > 0.0) || !std::isfinite(V0)) { cr.spa_used[c][pi] = -1; continue; }
                if (K <= KREML) {                 // small genes: REML-GLS test (below)
                    double tst, nu;
                    // sensitivity of the conditional mean to trait 1's tested component, and the
                    // sampling variance of its estimate (same trace formula as for trait 2)
                    prep_vth(t1);
                    std::vector<Eigen::VectorXd> msl(A);
                    const Eigen::VectorXd p1a = pcache[t1].solve(a1);
                    std::vector<Eigen::VectorXd> xbg;          // background covariates B_a Psi1^{-1} a1, a != c
                    for (int a2 = 0; a2 < A; ++a2) {
                        get_M(a2);
                        const Eigen::VectorXd mp = Mc[a2] * p1a;
                        msl[a2] = -(Psi12 * pcache[t1].solve(mp));
                        if (a2 != c && g[a2] != 0.0 && off[a2 + 1] > off[a2]) xbg.push_back(mp);
                    }
                    if (coher_reml_test(c, t1, t2, b, a1, r1, e1, e2, e12, w, Psi12, t2p.P0, Mc[c], pcache[t1],
                                        ec, yty, ycross[pi], msl, vth[t1], xbg, tst, nu)) {
                        cr.p_spa[c][pi] = t_two_sided_p(tst, nu); cr.spa_used[c][pi] = 7;
                    } else { cr.p_spa[c][pi] = NA_REAL; cr.spa_used[c][pi] = -1; }
                    continue;
                }
                // large genes: s_hat has many degrees of freedom; include it only when clearly
                // positive (> 2 null SD), with a Satterthwaite t for its uncertainty
                const double sd0 = std::sqrt(2.0 * t2p.tau0);
                double V = V0, nu = 1e300;
                if (s2 > 2.0 * sd0) {
                    const double var_s = 2.0 * (t2p.tau0 + t2p.tau1 * s2 + tau2c[c] * s2 * s2);
                    V = V0 + s2 * bm;
                    nu = 2.0 * V * V / (bm * bm * var_s);
                }
                cr.p_spa[c][pi] = t_two_sided_p((Q - hm) / std::sqrt(V), nu);
                cr.spa_used[c][pi] = 6;
            }
        }
        finish_se();
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

        // ---- exact binary SPA: the exact distribution of every SPARSE column's score ----
        // u_j = mu_j + delta * kappa_j * S_j, S_j = sum over carriers of dosage x case status.
        // A column is sparse when its expected number of case carriers is <= ec_max.
        // Sparse columns that SHARE carriers (MNP pieces, split multi-allelics, variants on
        // one haplotype) are not independent: they are grouped (connected components of
        // shared carriers) and each group gets ONE exact joint distribution. Within a group
        // the carriers fall into "atoms" (individuals with the same dosage in every column of
        // the group); the case counts of the atoms are independent Poisson-binomials and
        // every S_j is a fixed combination of them, so the joint null is enumerated exactly
        // (configurations with probability < 1e-14 are dropped). A group whose enumeration
        // exceeds 50,000 configurations falls back to independent columns.
        struct SpCol { std::vector<int> js; std::vector<double> U, lp, uobs; int obs; };   // U: npts x m
        std::vector<SpCol> spc; std::vector<int> dns; std::vector<double> dns_mean;
        if (use_exact) {
            const double dl = bx.delta[t], Ncase = bx.ncase[t];
            const std::vector<unsigned char>& bb = *bx.b[t];
            struct ColP { int j; double mu, kap; };
            std::vector<ColP> sp;
            for (int j = 0; j < K; ++j) {
                const std::vector<int>& ci = cd.car[j];
                const std::vector<unsigned char>& di = cd.dos[j];
                double ec = 0.0;
                for (size_t r = 0; r < ci.size(); ++r) ec += (double) bx.PiF(ci[r], t) * di[r];
                if (ci.empty() || ec > bx.ec_max) {
                    dns.push_back(j); dns_mean.push_back(Ulo(j, t) + dl * Upi(j, t)); continue;
                }
                const int cn = (int) ci.size();
                double scar = 0.0;
                for (int r = 0; r < cn; ++r) scar += cd.V(ci[r], j);
                const double xnc = (cd.colsum[j] - scar) / std::max(1, n - cn);   // non-carrier value
                double num = 0.0, den = 0.0;
                for (int r = 0; r < cn; ++r) { num += ((double) cd.V(ci[r], j) - xnc) * di[r]; den += (double) di[r] * di[r]; }
                ColP cp; cp.j = j; cp.kap = num / den; cp.mu = Ulo(j, t) + dl * xnc * Ncase;
                sp.push_back(cp);
            }
            // single column: Poisson-binomial of S_j
            auto add_single = [&](const ColP& cp) {
                const std::vector<int>& ci = cd.car[cp.j];
                const std::vector<unsigned char>& di = cd.dos[cp.j];
                const int cn = (int) ci.size();
                std::vector<double> pis(cn); int sobs = 0;
                for (int r = 0; r < cn; ++r) { pis[r] = bx.PiF(ci[r], t); if (bb[ci[r]]) sobs += di[r]; }
                SpCol sc; sc.js.assign(1, cp.j);
                pb_pmf(di, pis, sc.lp);
                if ((int) sc.lp.size() < sobs + 1) sc.lp.resize(sobs + 1, -1e300);   // keep the observation
                sc.U.resize(sc.lp.size());
                for (size_t q2 = 0; q2 < sc.U.size(); ++q2) sc.U[q2] = cp.mu + dl * cp.kap * (double) q2;
                sc.uobs.assign(1, cp.mu + dl * cp.kap * (double) sobs); sc.obs = sobs;
                spc.push_back(sc);
            };
            // connected components of shared carriers (union-find)
            const int ns = (int) sp.size();
            std::vector<int> par(ns);
            for (int a = 0; a < ns; ++a) par[a] = a;
            auto fnd = [&](int x) { while (par[x] != x) { par[x] = par[par[x]]; x = par[x]; } return x; };
            {
                // two columns are joined when their case counts are strongly correlated:
                // corr(S_j, S_k) = sum_shared d_j d_k / sqrt(sum d_j^2 sum d_k^2) >= 0.3 (identical
                // or nested carrier sets). Edges are added strongest first and a group is capped
                // at 16 columns, so chance overlaps cannot chain whole genes together.
                std::unordered_map<int, std::vector<std::pair<int, unsigned char> > > owner;
                std::vector<double> ss(ns, 0.0);
                for (int a = 0; a < ns; ++a) {
                    const int j = sp[a].j;
                    for (size_t r = 0; r < cd.car[j].size(); ++r) {
                        owner[cd.car[j][r]].push_back(std::make_pair(a, cd.dos[j][r]));
                        ss[a] += (double) cd.dos[j][r] * cd.dos[j][r];
                    }
                }
                std::map<std::pair<int, int>, double> shared;
                for (std::unordered_map<int, std::vector<std::pair<int, unsigned char> > >::iterator it = owner.begin();
                     it != owner.end(); ++it) {
                    const std::vector<std::pair<int, unsigned char> >& L2 = it->second;
                    for (size_t x = 0; x < L2.size(); ++x)
                        for (size_t y = x + 1; y < L2.size(); ++y)
                            shared[std::make_pair(std::min(L2[x].first, L2[y].first), std::max(L2[x].first, L2[y].first))]
                                += (double) L2[x].second * L2[y].second;
                }
                std::vector< std::pair<double, std::pair<int, int> > > edges;
                for (std::map<std::pair<int, int>, double>::iterator it = shared.begin(); it != shared.end(); ++it) {
                    const double rho = it->second / std::sqrt(ss[it->first.first] * ss[it->first.second]);
                    if (rho >= 0.3) edges.push_back(std::make_pair(-rho, it->first));
                }
                std::sort(edges.begin(), edges.end());
                std::vector<int> gsz(ns, 1);
                for (size_t e2 = 0; e2 < edges.size(); ++e2) {
                    const int r1 = fnd(edges[e2].second.first), r2 = fnd(edges[e2].second.second);
                    if (r1 == r2 || gsz[r1] + gsz[r2] > 16) continue;
                    par[r1] = r2; gsz[r2] += gsz[r1];
                }
            }
            std::map<int, std::vector<int> > grp;
            for (int a = 0; a < ns; ++a) grp[fnd(a)].push_back(a);
            for (std::map<int, std::vector<int> >::iterator gi = grp.begin(); gi != grp.end(); ++gi) {
                const std::vector<int>& mem = gi->second;
                if (mem.size() == 1) { add_single(sp[mem[0]]); continue; }
                const int m = (int) mem.size();
                // atoms: carriers with the same dosage vector over the group's columns
                std::map<int, std::vector<unsigned char> > dv;
                for (int k = 0; k < m; ++k) {
                    const int j = sp[mem[k]].j;
                    for (size_t r = 0; r < cd.car[j].size(); ++r) {
                        std::vector<unsigned char>& v = dv[cd.car[j][r]];
                        if (v.empty()) v.assign(m, 0);
                        v[k] = cd.dos[j][r];
                    }
                }
                std::map<std::vector<unsigned char>, int> pat;
                std::vector< std::vector<unsigned char> > apat; std::vector< std::vector<double> > api; std::vector<int> aobs;
                for (std::map<int, std::vector<unsigned char> >::iterator it = dv.begin(); it != dv.end(); ++it) {
                    std::map<std::vector<unsigned char>, int>::iterator pt = pat.find(it->second);
                    int a;
                    if (pt == pat.end()) { a = (int) apat.size(); pat[it->second] = a; apat.push_back(it->second); api.push_back(std::vector<double>()); aobs.push_back(0); }
                    else a = pt->second;
                    api[a].push_back(bx.PiF(it->first, t));
                    if (bb[it->first]) ++aobs[a];
                }
                const int na = (int) apat.size();
                std::vector< std::vector<double> > alp(na);
                std::vector<double> amax(na + 1, 0.0);
                for (int a = 0; a < na; ++a) {
                    std::vector<unsigned char> ones(api[a].size(), 1);
                    pb_pmf(ones, api[a], alp[a]);
                    if ((int) alp[a].size() < aobs[a] + 1) alp[a].resize(aobs[a] + 1, -1e300);
                }
                for (int a = na - 1; a >= 0; --a) {
                    double mx = -1e300; for (double v : alp[a]) mx = std::max(mx, v);
                    amax[a] = amax[a + 1] + mx;                      // best case for atoms a..na-1
                }
                const double LMIN = std::log(1e-14);
                const size_t CAP = 50000;
                SpCol sc; sc.js.resize(m); sc.obs = -1;
                for (int k = 0; k < m; ++k) sc.js[k] = sp[mem[k]].j;
                std::vector<int> cnt(na, 0); std::vector<double> Sk(m, 0.0);
                bool over = false;
                // depth-first enumeration of the atom case counts, pruned on probability
                std::function<void(int, double)> rec = [&](int a, double lp) {
                    if (over) return;
                    if (a == na) {
                        if (sc.lp.size() >= CAP) { over = true; return; }
                        bool isobs = true;
                        for (int b = 0; b < na; ++b) if (cnt[b] != aobs[b]) { isobs = false; break; }
                        if (isobs) sc.obs = (int) sc.lp.size();
                        sc.lp.push_back(lp);
                        for (int k = 0; k < m; ++k) sc.U.push_back(sp[mem[k]].mu + dl * sp[mem[k]].kap * Sk[k]);
                        return;
                    }
                    for (int c2 = 0; c2 < (int) alp[a].size(); ++c2) {
                        const double l2 = lp + alp[a][c2];
                        if (l2 + amax[a + 1] < LMIN) continue;
                        cnt[a] = c2;
                        for (int k = 0; k < m; ++k) Sk[k] += (double) apat[a][k] * c2;
                        rec(a + 1, l2);
                        for (int k = 0; k < m; ++k) Sk[k] -= (double) apat[a][k] * c2;
                        if (over) return;
                    }
                    cnt[a] = 0;
                };
                rec(0, 0.0);
                if (over || sc.lp.empty()) { for (int k = 0; k < m; ++k) add_single(sp[mem[k]]); continue; }
                std::vector<double> So(m, 0.0); double lpo = 0.0;
                for (int a = 0; a < na; ++a) {
                    lpo += alp[a][aobs[a]];
                    for (int k = 0; k < m; ++k) So[k] += (double) apat[a][k] * aobs[a];
                }
                sc.uobs.resize(m);
                for (int k = 0; k < m; ++k) sc.uobs[k] = sp[mem[k]].mu + dl * sp[mem[k]].kap * So[k];
                if (sc.obs < 0) {                                    // observation was pruned: keep it
                    sc.obs = (int) sc.lp.size(); sc.lp.push_back(std::max(lpo, -1e300));
                    for (int k = 0; k < m; ++k) sc.U.push_back(sc.uobs[k]);
                }
                spc.push_back(sc);
            }
        }
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

            // Every test goes to the saddlepoint (no Wald screen). The Wald p
            // from the closed-form variance is kept ONLY as the fallback for a
            // saddlepoint that fails (spa_used = 0).
            double p_wald = one_sided ? 0.5 * std::erfc((sigma[c] / std::sqrt(var0)) / std::sqrt(2.0))
                                      : std::erfc(std::abs(sigma[c] / std::sqrt(var0)) / std::sqrt(2.0));

            // ---- exact per-variant binary saddlepoint (spa_used = 5) -------------------
            if (use_exact) {
                const double Cc = Tinv(c, env) * yty[t];
                double qobs = Cc;
                std::vector<ExactTerm> terms; terms.reserve(spc.size());
                for (size_t r = 0; r < spc.size(); ++r) {
                    const SpCol& sc = spc[r];
                    const int m = (int) sc.js.size();
                    std::vector<double> ak(m); bool any = false;
                    for (int k = 0; k < m; ++k) { ak[k] = Tinv(c, comp_of[sc.js[k]]) * g[comp_of[sc.js[k]]]; if (ak[k] != 0.0) any = true; }
                    if (!any) continue;
                    ExactTerm e; e.lp = sc.lp; e.v.assign(sc.lp.size(), 0.0); e.obs = sc.obs;
                    for (size_t q2 = 0; q2 < e.v.size(); ++q2)
                        for (int k = 0; k < m; ++k) { const double u = sc.U[q2 * m + k]; e.v[q2] += ak[k] * u * u; }
                    for (int k = 0; k < m; ++k) qobs += ak[k] * sc.uobs[k] * sc.uobs[k];
                    terms.push_back(e);
                }
                std::vector<double> lamd, nc2d;
                bool ok_d = true;
                if (!dns.empty()) {
                    const int nd = (int) dns.size();
                    Eigen::MatrixXd LD(nd, K);
                    for (int r = 0; r < nd; ++r) LD.row(r) = L.row(dns[r]);
                    // Cov of the dense scores under the Bernoulli null: delta^2 X_D' diag(pi(1-pi)) X_D,
                    // approximated by vbar X_D' X_D (exact when pi is constant). NOT the Gaussian-null
                    // Smat: its residual variance is estimated from y (minus any LOCO score), which
                    // is not the Bernoulli variance when a polygenic score is subtracted.
                    Eigen::MatrixXd Sig = bx.vbar[t] * (LD * LD.transpose());
                    Sig = (0.5 * (Sig + Sig.transpose())).eval();
                    Eigen::VectorXd ad(nd), md(nd);
                    for (int r = 0; r < nd; ++r) {
                        const int j = dns[r];
                        ad[r] = Tinv(c, comp_of[j]) * g[comp_of[j]];
                        md[r] = dns_mean[r];
                        qobs += ad[r] * U(j, t) * U(j, t);
                    }
                    const double sscale = Sig.diagonal().mean();
                    Eigen::LLT<Eigen::MatrixXd> ch; double jit = 0.0; ok_d = false;
                    for (int it = 0; it < 6; ++it) {
                        Eigen::MatrixXd Sj = Sig; Sj.diagonal().array() += jit * sscale;
                        ch.compute(Sj);
                        if (ch.info() == Eigen::Success) { ok_d = true; break; }
                        jit = (jit == 0.0) ? 1e-10 : jit * 100.0;
                    }
                    if (ok_d) {
                        Eigen::MatrixXd Fd = ch.matrixL();
                        Eigen::MatrixXd Bd = Fd.transpose() * ad.asDiagonal() * Fd;
                        Bd = (0.5 * (Bd + Bd.transpose())).eval();
                        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Bd);
                        if (es.info() != Eigen::Success) ok_d = false;
                        else {
                            Eigen::VectorXd nu = es.eigenvectors().transpose() *
                                                 Fd.triangularView<Eigen::Lower>().solve(md);
                            for (int r = 0; r < nd; ++r) {
                                const double l = es.eigenvalues()[r];
                                if (std::abs(l) < 1e-300) continue;
                                lamd.push_back(l); nc2d.push_back(nu[r] * nu[r]);
                            }
                        }
                    }
                }
                if (ok_d && (!terms.empty() || !lamd.empty())) {
                    if (g_exact_capture && !g_exact_capture->done && c == 0) {
                        g_exact_capture->T = terms; g_exact_capture->lam = lamd; g_exact_capture->nc2 = nc2d;
                        g_exact_capture->C = Cc; g_exact_capture->q = qobs; g_exact_capture->done = true;
                    }
                    ExactSpaResult er = exact_spa_solve(qobs, terms, lamd, nc2d, Cc);
                    // saddlepoint failed, or p small: lattice convolution (see exact_grid_solve)
                    if (!terms.empty() && !er.exact &&
                        (!er.converged || (one_sided ? er.p_up : std::min(er.p_up, er.p_lo)) < bx.grid_p)) {
                        ExactSpaResult eg = exact_grid_solve(qobs, terms, lamd, nc2d, Cc);
                        if (eg.converged) er = eg;
                    }
                    if (er.converged) {
                        cr.p_spa[c][t] = one_sided ? er.p_up : std::min(1.0, 2.0 * std::min(er.p_up, er.p_lo));
                        cr.spa_used[c][t] = 5;
                        continue;
                    }
                }
            }
            // binary phenotypes use the exact computation only: if it is not possible, the
            // p value is reported as NA (spa_used = -1) rather than from an approximation
            if (binary) { cr.p_spa[c][t] = NA_REAL; cr.spa_used[c][t] = -1; continue; }


            // ---- quantitative: build Asym and solve the eigenproblem ---------
            Eigen::LLT<Eigen::MatrixXd> lltS(Smat);
            if (lltS.info() != Eigen::Success) { cr.p_spa[c][t] = p_wald; cr.spa_used[c][t] = 0; continue; }
            Eigen::MatrixXd Cf = lltS.matrixL();
            Eigen::MatrixXd Asym = Cf.transpose() * Amat * Cf;
            Asym = (0.5 * (Asym + Asym.transpose())).eval();
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ses(Asym, Eigen::EigenvaluesOnly);
            if (ses.info() != Eigen::Success) { cr.p_spa[c][t] = p_wald; cr.spa_used[c][t] = 0; continue; }
            const Eigen::VectorXd& ev = ses.eigenvalues();
            std::vector<double> eig(ev.data(), ev.data() + K);
            QuadSpaResult qr = quad_spa_solve(sigma[c], eig, eig_rep, n_rep, 100, 1e-8,
                                              0.0, 0.0, 0.0, 0.0, one_sided ? 1 : 0);
            if (qr.converged) { cr.p_spa[c][t] = qr.p; cr.spa_used[c][t] = 1; }
            else              { cr.p_spa[c][t] = p_wald; cr.spa_used[c][t] = 0; }
        }
    }
    finish_se();
}


// ===========================================================================
// Driver
// ===========================================================================
struct ChunkParams {
    int chunk_size; long window_bp; bool do_windows; int max_window_chunks;
    int max_gene_snps, min_gene_snps;
    bool spa; bool binary;
    double acatv_mac;          // variants with MAC below this are collapsed into a burden
    bool h2_raw;               // with loco_prs: vg / se / h2 from y itself, p-values from y - LOCO
    bool one_sided;            // univariate tests: upper-tail p (variance components >= 0)
    int  cov_df;
    std::string out_file; int batch_size, n_threads;
    std::string chr;
    bool coher;
    std::vector< std::pair<int,int> > pairs;
    std::vector<std::string> pair_names;
    // conditioning of gene tests on the gene's lead common variants (see cond_common_apply)
    bool cond_common = false; double cond_p = 1e-4; int cond_max = 10; double cond_maf = 0.01;
    bool use_acatv = false;    // ACAT-V (marginal single-variant tests): off with coher or cond_common
};

// One unit of work: its own range [a,b) plus two flank ranges. `win` is the
// index of the window a pass-1 chunk belongs to (or of the window itself).
struct Job { size_t ci; int a, b, fL0, fL1, fR0, fR1; int win; bool gene = false; };

// ---------------------------------------------------------------------------
// ACAT-V for one category of one unit (gene block / chunk), at one alpha.
//
// Each variant contributes its KVIK p-value p_j = erfc(|Z_j|/sqrt2) to a
// Cauchy combination with weight
//     w_j = [2 f_j (1 - f_j)]^(1 + alpha)
// -- the variant's expected share of heritability under the alpha model, the
// same quantity that weights it in the variance-component test. (The ACAT-V
// default, Beta(f;1,25)^2 f(1-f), is ~ 2f(1-f) for rare variants, i.e. alpha = 0.)
//
// As in ACAT-V, variants with minor allele count below `macthr` are not
// combined individually but collapsed into one burden statistic,
//     Z_b = sum_j s_j Z_j / sqrt(sum_j s_j^2),   s_j = sqrt(w_j),
// with Z_j oriented to the minor allele. This uses the KVIK statistics as if
// the collapsed variants were uncorrelated, which holds closely for ultra-rare
// variants (they almost never share carriers). The burden enters the Cauchy
// combination as one variant, with the average weight of its members.
// ---------------------------------------------------------------------------
static double acatv_p(const ChunkDataA& cd, size_t qcat,
                      const std::vector<float>& kz, const std::vector<float>& kmac,
                      double alpha, double macthr) {
    std::vector<double> pv, pw;
    double bz = 0.0, bs2 = 0.0; int bn = 0;
    if (qcat >= cd.cat_obim.size()) return NA_REAL;
    const std::vector<int>& vb = cd.cat_obim[qcat];
    const std::vector<float>& vf = cd.cat_omaf[qcat];
    for (size_t j = 0; j < vb.size(); ++j) {
        const int b = vb[j];
        if (b < 0) continue;
        const double z = kz[b];
        if (std::isnan(z)) continue;
        const double f = vf[j];
        if (!(f > 0.0 && f < 1.0)) continue;
        const double s2 = (alpha == -1.0) ? 1.0 : std::pow(2.0 * f * (1.0 - f), 1.0 + alpha);
        if (!(s2 > 0.0) || !std::isfinite(s2)) continue;
        if (kmac[b] >= macthr) {
            pv.push_back(std::erfc(std::abs(z) / std::sqrt(2.0)));
            pw.push_back(s2);
        } else {
            bz += std::sqrt(s2) * z; bs2 += s2; ++bn;
        }
    }
    if (bn > 0 && bs2 > 0.0) {
        const double zb = bz / std::sqrt(bs2);
        pv.push_back(std::erfc(std::abs(zb) / std::sqrt(2.0)));
        pw.push_back(bs2 / bn);
    }
    if (pv.empty()) return NA_REAL;
    return acat_combine(pv, pw);
}

// The worker takes a VECTOR of phenotype variants (distinct prs_mask rows).
// Each job is read from disk and its geometry built ONCE; the variants are then
// swept over that geometry. Pass 1 passes a single variant, so it is unaffected.
struct Worker : public RcppParallel::Worker {
    const ChunkContext& ctx; const std::vector<Job>& jobs; size_t job0;
    const std::vector<PhenoStats>& S; const ChunkParams& pr; bool spa;
    std::vector< std::vector<ChunkResultA> >& out;      // [variant][slot]
    std::vector<char>& ok;
    Worker(const ChunkContext& ctx, const std::vector<Job>& jobs, size_t job0,
           const std::vector<PhenoStats>& S, const ChunkParams& pr, bool spa,
           std::vector< std::vector<ChunkResultA> >& out, std::vector<char>& ok)
        : ctx(ctx), jobs(jobs), job0(job0), S(S), pr(pr), spa(spa), out(out), ok(ok) {}
    // Column weight exactly as apply_alpha() applies it, so the ratio between two
    // alphas reproduces what reading the data at the second alpha would give.
    static double alpha_w(double f, double a) {
        if (a == -1.0) return 1.0;
        return (f > 0.0 && f < 1.0) ? std::pow(2.0 * f * (1.0 - f), (1.0 + a) / 2.0) : 1.0;
    }
    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t w = begin; w < end; ++w) {
            const Job& j = jobs[job0 + w];
            ChunkDataA cd;
            if (!make_chunk(ctx, j.ci, j.a, j.b, j.fL0, j.fL1, j.fR0, j.fR1, cd)) { ok[w] = 0; continue; }
            // cond_common: gene blocks are tested conditional on their lead common variants
            Eigen::MatrixXd Qc; std::string cnames; int rc = 0;
            std::vector<PhenoStats> Sloc;
            if (pr.cond_common && j.gene) {
                CondOpt co; co.p = pr.cond_p; co.maf = pr.cond_maf; co.max = pr.cond_max;
                rc = cond_common_apply(ctx, co, cd, Qc, cnames, S.empty() ? 0 : &S[0].Y);
                if (rc < 0) { ok[w] = 0; continue; }
                if (rc > 0) {
                    Sloc.resize(S.size());
                    for (size_t v = 0; v < S.size(); ++v)
                        cond_pheno_stats(ctx, pr.binary, pr.coher, pr.pairs, S[v], Qc, Sloc[v]);
                }
            }
            const std::vector<PhenoStats>& SS = (rc > 0) ? Sloc : S;
            const int cdf = pr.cov_df + std::max(rc, 0);
            // Multi-alpha applies to ASSOCIATION passes only (spa = true). The
            // heritability pass runs with spa = false and stays single-alpha.
            const bool multi = spa && ctx.alphas.size() > 1;
            // With h2_raw, S = { y - LOCO (p-values), y (heritability) }. The extra
            // alphas and ACAT-V only concern p-values, so they run on variant 0 only;
            // variant 1 is tested once, at the primary alpha, which is all its
            // reported se_vg needs.
            const size_t NVP = (spa && pr.h2_raw) ? std::min<size_t>(1, SS.size()) : SS.size();
            ChunkGeom gm; Eigen::MatrixXd G0;
            build_geom(cd, cdf, spa, pr.binary, gm, 0, multi ? &G0 : 0);   // may leave gm.ok false
            for (size_t v = 0; v < SS.size(); ++v)
                test_chunk_annot(cd, gm, SS[v].Y, SS[v].Yf, SS[v].Vp, SS[v].yty,
                                 pr.binary, spa, pr.coher, pr.pairs,
                                 SS[v].ycross, SS[v].cb, SS[v].bx, pr.one_sided, out[v][w]);
            if (multi) {
                // ---- ACAT over alpha ------------------------------------------
                // V was read at alphas[0]. Another alpha only rescales columns:
                //     V_a = V_0 D,   d_j = w(f_j, a) / w(f_j, a0)
                // so  G_a = D G_0 D  and  ||v_j||^2 scales by d_j^2. No second
                // read, no second n K^2 Gram -- the extra cost per alpha is one
                // copy of V for the n x K GEMM V'Y and the per-test SPA.
                const size_t NAL = ctx.alphas.size();
                const double a0 = ctx.alphas[0];
                for (size_t v = 0; v < NVP; ++v) {
                    ChunkResultA& R = out[v][w];
                    R.p_alpha.assign(R.p_spa.size(), std::vector< std::vector<double> >());
                    for (size_t q = 0; q < R.p_spa.size(); ++q) {
                        R.p_alpha[q].assign(R.p_spa[q].size(), std::vector<double>(NAL, NA_REAL));
                        for (size_t t = 0; t < R.p_spa[q].size(); ++t) R.p_alpha[q][t][0] = R.p_spa[q][t];
                    }
                }
                if (gm.ok) {
                    for (size_t k = 1; k < NAL; ++k) {
                        const double ak = ctx.alphas[k];
                        Eigen::VectorXd d(cd.K);
                        for (int c = 0; c < cd.K; ++c)
                            d[c] = alpha_w(cd.colmaf[c], ak) / alpha_w(cd.colmaf[c], a0);
                        ChunkDataA ck = cd;
                        for (int c = 0; c < cd.K; ++c) {
                            ck.V.col(c) *= (float) d[c];
                            ck.cj[c] *= d[c] * d[c];
                        }
                        Eigen::MatrixXd Gk = d.asDiagonal() * G0 * d.asDiagonal();
                        ChunkGeom gk;
                        build_geom(ck, cdf, spa, pr.binary, gk, &Gk);
                        for (size_t v = 0; v < NVP; ++v) {
                            ChunkResultA tmp;
                            test_chunk_annot(ck, gk, SS[v].Y, SS[v].Yf, SS[v].Vp, SS[v].yty,
                                             pr.binary, spa, pr.coher, pr.pairs,
                                             SS[v].ycross, SS[v].cb, SS[v].bx, pr.one_sided, tmp);
                            ChunkResultA& R = out[v][w];
                            for (size_t q = 0; q < R.p_alpha.size() && q < tmp.p_spa.size(); ++q)
                                for (size_t t = 0; t < R.p_alpha[q].size() && t < tmp.p_spa[q].size(); ++t)
                                    R.p_alpha[q][t][k] = tmp.p_spa[q][t];
                        }
                    }
                }
                // combine; vg / se_vg / h2 / spa_used stay at the primary alpha
                for (size_t v = 0; v < NVP; ++v) {
                    ChunkResultA& R = out[v][w];
                    for (size_t q = 0; q < R.p_alpha.size(); ++q)
                        for (size_t t = 0; t < R.p_alpha[q].size(); ++t) {
                            std::vector<double> pv, pw;
                            for (size_t k = 0; k < NAL; ++k)
                                if (!std::isnan(R.p_alpha[q][t][k])) { pv.push_back(R.p_alpha[q][t][k]); pw.push_back(1.0); }
                            R.p_spa[q][t] = pv.empty() ? NA_REAL : acat_combine(pv, pw);
                        }
                }
            }
            // ---- ACAT-V from LDAK-KVIK single-variant statistics -------------
            // Association passes only, univariate only. For every category and
            // trait: one ACAT-V p-value per alpha, then the omnibus
            //     p_spa = ACAT( VC p at each alpha , ACAT-V p at each alpha )
            // with equal weights (2 x n_alpha p-values). p_vc and p_acv keep the
            // two families separately (each ACAT-combined over alpha).
            if (spa && pr.use_acatv && !ctx.kz.empty()) {
                const size_t NAL = ctx.alphas.size();
                for (size_t v = 0; v < NVP; ++v) {
                    ChunkResultA& R = out[v][w];
                    const size_t A = R.cat_m.size();
                    const size_t NT = R.p_spa.empty() ? 0 : R.p_spa[0].size();
                    R.p_acatv.assign(A, std::vector< std::vector<double> >(NT, std::vector<double>(NAL, NA_REAL)));
                    R.p_vc.assign(A, std::vector<double>(NT, NA_REAL));
                    R.p_acv.assign(A, std::vector<double>(NT, NA_REAL));
                    for (size_t q = 0; q < A; ++q) {
                        for (size_t t = 0; t < NT && t < ctx.kz.size(); ++t) {
                            std::vector<double> av, aw, ov, ow;
                            for (size_t k = 0; k < NAL; ++k) {
                                const double pa = acatv_p(cd, q, ctx.kz[t], ctx.kmac[t],
                                                          ctx.alphas[k], pr.acatv_mac);
                                R.p_acatv[q][t][k] = pa;
                                if (!std::isnan(pa)) { av.push_back(pa); aw.push_back(1.0); }
                            }
                            R.p_vc[q][t]  = R.p_spa[q][t];     // already the ACAT over alphas
                            R.p_acv[q][t] = av.empty() ? NA_REAL : acat_combine(av, aw);
                            if (!R.p_alpha.empty() && q < R.p_alpha.size() && t < R.p_alpha[q].size()) {
                                for (size_t k = 0; k < NAL; ++k)
                                    if (!std::isnan(R.p_alpha[q][t][k])) { ov.push_back(R.p_alpha[q][t][k]); ow.push_back(1.0); }
                            } else if (!std::isnan(R.p_spa[q][t])) { ov.push_back(R.p_spa[q][t]); ow.push_back(1.0); }
                            for (size_t k = 0; k < av.size(); ++k) { ov.push_back(av[k]); ow.push_back(1.0); }
                            R.p_spa[q][t] = ov.empty() ? NA_REAL : acat_combine(ov, ow);
                        }
                    }
                }
            }
            for (size_t v = 0; v < out.size(); ++v) { out[v][w].n_cond = rc; out[v][w].cond_snps = cnames; }
            ok[w] = 1;
        }
    }
};

static void run_jobs(const ChunkContext& ctx, const ChunkParams& pr,
                     const std::vector<PhenoStats>& S,
                     bool spa, const std::vector<Job>& jobs, size_t j0, size_t j1,
                     std::vector< std::vector<ChunkResultA> >& res, std::vector<char>& ok,
                     const char* tag, long& done_total, long total, int batch,
                     const std::chrono::steady_clock::time_point& t0) {
    typedef std::chrono::steady_clock clk;
    const size_t NV = S.size();
    size_t idx = j0;
    while (idx < j1) {
        const size_t nb = std::min((size_t) std::max(1, batch), j1 - idx);
        std::vector< std::vector<ChunkResultA> > tmp(NV, std::vector<ChunkResultA>(nb));
        std::vector<char> tok(nb, 0);
        Worker wk(ctx, jobs, idx, S, pr, spa, tmp, tok);
        RcppParallel::parallelFor(0, nb, wk);
        for (size_t v = 0; v < NV; ++v)
            for (size_t b = 0; b < nb; ++b) res[v][idx + b] = tmp[v][b];
        for (size_t b = 0; b < nb; ++b) ok[idx + b] = tok[b];
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
        int wa = lo, wn = 0;                     // window start, chunks in it so far
        for (int a = lo; a <= hi; ) {
            const int b = std::min(hi + 1, a + pr.chunk_size);
            Job j; j.ci = ci; j.a = a; j.b = b; j.fL0 = j.fL1 = j.fR0 = j.fR1 = 0;
            j.win = (int) wj.size() + (int) wr.size();
            cj.push_back(j); ++wn;
            const bool last = (b > hi);
            // Windows close on bp span OR on the chunk cap, whichever comes
            // first. Without the cap a gene-dense megabase can reach tens of
            // thousands of columns, and both the Gram (nK^2/2) and the memory
            // (n*K floats for V, K^2 for the Gram) blow up with it.
            if (last || wn >= pr.max_window_chunks ||
                ctx.bim.bp[b - 1] - ctx.bim.bp[wa] + 1 >= pr.window_bp) {
                wr.push_back(std::make_pair(wa, b)); wa = b; wn = 0;
            }
            a = b;
        }
        // A trailing remainder shorter than half a window is absorbed into its
        // predecessor (as the old build_cells did), so no window-level estimate
        // rests on a sliver. The chunks it held are re-pointed accordingly.
        // Skipped when the merge would push the predecessor past the cap.
        if (wr.size() >= 2) {
            const std::pair<int,int>& lastw = wr.back();
            const int merged = ((lastw.second - wr[wr.size() - 2].first)
                                + pr.chunk_size - 1) / pr.chunk_size;
            if (ctx.bim.bp[lastw.second - 1] - ctx.bim.bp[lastw.first] + 1 < pr.window_bp / 2
                && merged <= pr.max_window_chunks) {
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
        int Kmax = 0; double Kmean = 0.0;
        for (size_t w = 0; w < wj.size(); ++w) {
            const int Kw = (wj[w].b - wj[w].a) + (wj[w].fL1 - wj[w].fL0) + (wj[w].fR1 - wj[w].fR0);
            Kmax = std::max(Kmax, Kw); Kmean += Kw;
        }
        if (!wj.empty()) Kmean /= (double) wj.size();
        Rcout << "Pass 2: " << wj.size() << " windows of ~" << pr.window_bp
              << " bp (at most " << pr.max_window_chunks << " chunks) with neighbouring "
              << "windows as flanks, phenotype = y - full PRS per prs_mask.\n";
        Rcout << "  mean K ~ " << (int) Kmean << ", worst case " << Kmax
              << "; cost per window scales as K^2.\n";
        report_plan(ctx, pr, Kmax, "pass 2");
    }

    // ---- GENE JOBS ---------------------------------------------------------
    // One job per gene, or several when a gene exceeds max_gene_snps. Genes are
    // tested marginally with NO flanks, matching SAIGE-GENE. They may overlap
    // and nest; nothing downstream assumes they partition anything, which is
    // exactly why gene results and the agnostic screen stay in separate files.
    struct GeneRec { std::string name, chr; long start, end; size_t j0, j1; };
    std::vector<Job> gj; std::vector<GeneRec> grec;
    std::vector<size_t> chr_g0, chr_g1;
    const bool has_genes = !ctx.gene_name.empty();
    long n_gene_drop_small = 0, n_gene_drop_chr = 0;
    std::map<std::string, long> gene_chr_unmatched;
    if (has_genes) {
        std::map<std::string, size_t> chr_of;
        for (size_t ci = 0; ci < ctx.chr_order.size(); ++ci)
            chr_of[strip_chr(ctx.chr_order[ci])] = ci;
        // Resolve each gene's chromosome ONCE, so the diagnostic below is exact
        // regardless of which chromosomes are being analysed.
        std::vector<int> gci(ctx.gene_name.size(), -1);
        for (size_t g = 0; g < ctx.gene_name.size(); ++g) {
            std::map<std::string, size_t>::const_iterator it = chr_of.find(strip_chr(ctx.gene_chr[g]));
            if (it == chr_of.end()) { ++gene_chr_unmatched[ctx.gene_chr[g]]; ++n_gene_drop_chr; }
            else gci[g] = (int) it->second;
        }
        for (size_t ci = 0; ci < ctx.chr_order.size(); ++ci) {
            if (!pr.chr.empty() && strip_chr(ctx.chr_order[ci]) != strip_chr(pr.chr)) continue;
            chr_g0.push_back(gj.size());
            for (size_t g = 0; g < ctx.gene_name.size(); ++g) {
                if (gci[g] != (int) ci) continue;
                const int lo = ctx.chr_lo[ci], hi = ctx.chr_hi[ci];
                const int a = lower_index(ctx.bim.bp, lo, hi, ctx.gene_start[g]);
                const int b = lower_index(ctx.bim.bp, lo, hi, ctx.gene_end[g] + 1);
                if (b - a < pr.min_gene_snps) { ++n_gene_drop_small; continue; }
                GeneRec r; r.name = ctx.gene_name[g]; r.chr = ctx.chr_order[ci];
                r.start = ctx.gene_start[g]; r.end = ctx.gene_end[g]; r.j0 = gj.size();
                for (int s = a; s < b; ) {
                    const int e = std::min(b, s + pr.max_gene_snps);
                    Job j; j.ci = ci; j.a = s; j.b = e;
                    j.fL0 = j.fL1 = j.fR0 = j.fR1 = 0;
                    j.win = (int) grec.size(); j.gene = true;
                    gj.push_back(j); s = e;
                }
                r.j1 = gj.size(); grec.push_back(r);
            }
            chr_g1.push_back(gj.size());
        }
        Rcout << "Genes: " << grec.size() << " with >= " << pr.min_gene_snps << " SNPs -> "
              << gj.size() << " blocks (split at " << pr.max_gene_snps << " SNPs), no flanks, "
              << "phenotype = y - LOCO PRS.\n";
        if (n_gene_drop_small) Rcout << "  " << n_gene_drop_small << " genes skipped: fewer than "
                                     << pr.min_gene_snps << " SNPs in the genotype file.\n";
        if (n_gene_drop_chr) {
            Rcout << "  " << n_gene_drop_chr << " genes skipped: chromosome not in the genotype file.\n";
            // Quote the labels: that is what makes stray whitespace visible.
            Rcout << "    unmatched gene chromosomes:";
            int shown = 0;
            for (std::map<std::string, long>::const_iterator it = gene_chr_unmatched.begin();
                 it != gene_chr_unmatched.end(); ++it) {
                if (shown++ == 8) { Rcout << " ..."; break; }
                Rcout << " '" << it->first << "'(" << it->second << ")";
            }
            Rcout << "\n    .bim chromosomes:";
            for (size_t ci = 0; ci < ctx.chr_order.size() && ci < 12; ++ci)
                Rcout << " '" << ctx.chr_order[ci] << "'";
            if (ctx.chr_order.size() > 12) Rcout << " ...";
            Rcout << "\n";
        }
        if (grec.empty()) stop("No gene overlapped the genotype file -- check the chromosome labels above, and the genome build of the gene coordinates");
    }
    // A genes-only run (genes supplied, do_windows = FALSE) skips the agnostic
    // chunk pass entirely: the user wants gene p-values, not chunk p-values.
    // With do_windows = TRUE the chunks are still needed, because the window
    // correction is defined as (window - sum of its chunks).
    const bool need_chunks = pr.do_windows || !has_genes;

    // ---- the per-chromosome LOCO offset, shared by the chunk and gene passes -
    const size_t G = ctx.prs_mask.empty() ? 0 : ctx.prs_mask[0].size();
    auto loco_stats = [&](size_t ci, std::vector<PhenoStats>& S) {
        Eigen::MatrixXd off = Eigen::MatrixXd::Zero(n, P);
        // ctx.prs is stored trait-major: prs[t * G + g]
        for (int t = 0; t < P; ++t)
            for (size_t g = 0; g < G; ++g) {
                const PrsSource& src = ctx.prs[(size_t) t * G + g];
                int col = -1;
                for (size_t c = 0; c < src.chr_index.size(); ++c) if (src.chr_index[c] == (int) ci) { col = (int) c; break; }
                if (col >= 0) off.col(t) += src.loco.col(col);
                else          off.col(t) += src.full;     // chromosome absent from the fit
            }
        S.resize(pr.h2_raw ? 2 : 1);
        build_pheno_stats(ctx, off, pr.binary, pr.coher, pr.pairs, S[0]);
        if (pr.h2_raw)                                    // the phenotype itself, no PRS
            build_pheno_stats(ctx, Eigen::MatrixXd(), pr.binary, pr.coher, pr.pairs, S[1]);
    };
    // HERITABILITY FROM y, P-VALUES FROM y - LOCO.
    // The LOCO PRS is fitted in-sample: y_hat = H y with H the ridge hat matrix of
    // the other chromosomes. Subtracting it removes a fraction ~ tr(H)/n of EVERY
    // component of y, including the tested region's own signal, so moment
    // estimates of vg shrink by ~ (1 - tr(H)/n)^2 (about 15-20% in simulation).
    // Tests stay valid and gain power (less residual noise), so p-values keep using
    // y - LOCO; vg / h2 / se_vg (and vg_flank, vg_env, vg_t1, vg_t2) are taken from
    // the same unit analysed on y itself.
    auto merge_raw = [&](std::vector< std::vector<ChunkResultA> >& rv, const std::vector<char>& okv) {
        if (!pr.h2_raw || rv.size() < 2) return;
        for (size_t w = 0; w < rv[0].size(); ++w) {
            if (!okv[w]) continue;
            ChunkResultA& R = rv[0][w]; const ChunkResultA& Q = rv[1][w];
            R.vg = Q.vg; R.h2 = Q.h2; R.se_vg = Q.se_vg;
            R.vg_flank = Q.vg_flank; R.vg_env = Q.vg_env;
            R.vg_t1 = Q.vg_t1; R.vg_t2 = Q.vg_t2;
        }
        std::vector<ChunkResultA>().swap(rv[1]);          // free the raw copies
    };
    const size_t NV1 = pr.h2_raw ? 2 : 1;

    // ---- PASS 1: association, one phenotype variant per chromosome ---------
    std::vector< std::vector<ChunkResultA> > r1v(NV1, std::vector<ChunkResultA>(need_chunks ? cj.size() : 0));
    std::vector<ChunkResultA>& r1 = r1v[0];
    std::vector<char> ok1(need_chunks ? cj.size() : 0, 0);
    if (need_chunks) {
        long done = 0;
        const std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
        for (size_t k = 0; k < chr_c0.size(); ++k) {
            if (chr_c1[k] == chr_c0[k]) continue;
            std::vector<PhenoStats> S; loco_stats(cj[chr_c0[k]].ci, S);
            run_jobs(ctx, pr, S, pr.spa, cj, chr_c0[k], chr_c1[k], r1v, ok1, "pass1",
                     done, (long) cj.size(), pr.batch_size, t1);
        }
        merge_raw(r1v, ok1);
        Rcout << "\n";
    } else {
        Rcout << "Skipping the agnostic chunk pass (genes supplied with do_windows = FALSE).\n";
    }

    // ---- PASS 1G: the same, on genes ---------------------------------------
    std::vector< std::vector<ChunkResultA> > rgv(NV1, std::vector<ChunkResultA>(gj.size()));
    std::vector<ChunkResultA>& rg = rgv[0];
    std::vector<char> okg(gj.size(), 0);
    if (has_genes && !gj.empty()) {
        long doneg = 0;
        const std::chrono::steady_clock::time_point tg = std::chrono::steady_clock::now();
        for (size_t k = 0; k < chr_g0.size(); ++k) {
            if (chr_g1[k] == chr_g0[k]) continue;
            std::vector<PhenoStats> S; loco_stats(gj[chr_g0[k]].ci, S);
            run_jobs(ctx, pr, S, pr.spa, gj, chr_g0[k], chr_g1[k], rgv, okg, "genes",
                     doneg, (long) gj.size(), pr.batch_size, tg);
        }
        merge_raw(rgv, okg);
        Rcout << "\n";
    }

    // ---- PASS 2: windows, one phenotype variant per distinct prs_mask row ---
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
        // All variants are built up front and swept inside the worker, so each
        // window is read and its Gram formed once rather than once per variant.
        std::vector<PhenoStats> S(variants.size());
        for (size_t v = 0; v < variants.size(); ++v) {
            Eigen::MatrixXd off = Eigen::MatrixXd::Zero(n, P);
            for (int t = 0; t < P; ++t)
                for (size_t g = 0; g < G; ++g)
                    if (variants[v][g]) off.col(t) += ctx.prs[(size_t) t * G + g].full;
            build_pheno_stats(ctx, off, false, pr.coher, pr.pairs, S[v]);
        }
        std::vector<char> okw(wj.size(), 0);
        long done2 = 0; const long tot2 = (long) wj.size();
        const std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
        const int b2 = std::max(1, std::min(pr.batch_size, 4));   // windows are slow: report often
        run_jobs(ctx, pr, S, false, wj, 0, wj.size(), r2, okw, "pass2", done2, tot2, b2, t2);
        for (size_t v = 0; v < variants.size(); ++v) ok2[v] = okw;
        Rcout << "\n";
    }

    // ---- ADJUSTMENT: distribute (window - sum of chunks) over the chunks -----
    // adj[chunk][k][t] for active category k of that chunk; NA if no window value
    // Under coher the two univariate components are adjusted by the SAME rule, so
    // a local genetic correlation can be formed with a calibrated numerator AND
    // calibrated denominators rather than mixing the two scales.
    // Sized from ok1, NOT from cj: in a genes-only run (need_chunks == false)
    // ok1 and r1 are deliberately empty, and indexing them by a chunk index
    // dereferences a null data pointer -- a segfault at address 0x0 right after
    // the gene pass finishes, which is exactly what that looks like.
    const size_t NCHUNK = ok1.size();
    std::vector< std::vector< std::vector<double> > > adj(NCHUNK), adj1(NCHUNK), adj2(NCHUNK);
    for (size_t k = 0; k < NCHUNK; ++k)
        if (ok1[k]) {
            const size_t nq = r1[k].cat_m.size();
            adj[k].assign(nq, std::vector<double>(NOUT, NA_REAL));
            if (pr.coher) {
                adj1[k].assign(nq, std::vector<double>(NOUT, NA_REAL));
                adj2[k].assign(nq, std::vector<double>(NOUT, NA_REAL));
            }
        }
    // per-unit-weight correction for each window x category x trait, reused by
    // the gene writer; NA where the window produced no estimate
    std::vector< std::vector< std::vector<double> > > delta(
        wj.size(), std::vector< std::vector<double> >(ctx.n_cat, std::vector<double>(NOUT, NA_REAL)));
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
                std::vector<double> ssum(NOUT, 0.0), ssum1(NOUT, 0.0), ssum2(NOUT, 0.0);
                std::vector< std::pair<size_t,int> > hits;
                for (size_t m = 0; m < members[w].size(); ++m) {
                    const size_t k = members[w][m]; const ChunkResultA& R = r1[k];
                    for (size_t q = 0; q < R.cat_id.size(); ++q) if (R.cat_id[q] == c) {
                        hits.push_back(std::make_pair(k, (int) q)); wsum += R.cat_w[q];
                        for (int t = 0; t < NOUT; ++t) {
                            if (!std::isnan(R.vg[q][t])) ssum[t] += R.vg[q][t];
                            if (pr.coher) {
                                if (!std::isnan(R.vg_t1[q][t])) ssum1[t] += R.vg_t1[q][t];
                                if (!std::isnan(R.vg_t2[q][t])) ssum2[t] += R.vg_t2[q][t];
                            }
                        }
                    }
                }
                if (hits.empty() || !(wsum > 0.0)) continue;
                // delta_w: the correction PER UNIT of alpha-model weight. The
                // chunks below get delta_w * cat_w; a gene inside this window
                // gets the same rule applied to its own cat_w (see the gene
                // writer), so both are calibrated by one definition.
                for (int t = 0; t < NOUT; ++t) {
                    const double Wv = W.vg[kw][t];
                    if (!std::isnan(Wv)) delta[w][c][t] = (Wv - ssum[t]) / wsum;
                }
                for (size_t h = 0; h < hits.size(); ++h) {
                    const size_t k = hits[h].first; const int q = hits[h].second;
                    const double share = r1[k].cat_w[q] / wsum;
                    for (int t = 0; t < NOUT; ++t) {
                        const double Wv = W.vg[kw][t];
                        if (!std::isnan(Wv) && !std::isnan(r1[k].vg[q][t]))
                            adj[k][q][t] = r1[k].vg[q][t] + (Wv - ssum[t]) * share;
                        if (!pr.coher) continue;
                        const double W1 = W.vg_t1[kw][t], W2 = W.vg_t2[kw][t];
                        if (!std::isnan(W1) && !std::isnan(r1[k].vg_t1[q][t]))
                            adj1[k][q][t] = r1[k].vg_t1[q][t] + (W1 - ssum1[t]) * share;
                        if (!std::isnan(W2) && !std::isnan(r1[k].vg_t2[q][t]))
                            adj2[k][q][t] = r1[k].vg_t2[q][t] + (W2 - ssum2[t]) * share;
                    }
                }
            }
        }
    }

    // ---- WRITE ---------------------------------------------------------------
    auto wr = [](std::ofstream& f, double v) { if (std::isnan(v)) f << "NA"; else f << v; };
    // multi-alpha: one extra column per alpha, p_alpha_<value>; p_spa = their ACAT
    const bool multi_a = pr.spa && ctx.alphas.size() > 1;
    std::string alpha_hdr;
    if (multi_a)
        for (size_t k = 0; k < ctx.alphas.size(); ++k) {
            std::ostringstream o; o << "\tp_alpha_" << ctx.alphas[k]; alpha_hdr += o.str();
        }
    // ACAT-V columns (KVIK summaries supplied): p_vc, p_acatv, p_acatv_<alpha>
    const bool acv = pr.use_acatv && !ctx.kz.empty();
    std::string acv_hdr;
    if (acv) {
        acv_hdr = "\tp_vc\tp_acatv";
        for (size_t k = 0; k < ctx.alphas.size(); ++k) {
            std::ostringstream o; o << "\tp_acatv_" << ctx.alphas[k]; acv_hdr += o.str();
        }
    }
    auto wr_acv = [&](std::ofstream& f, const ChunkResultA& R, size_t q, int t) {
        if (!acv) return;
        const bool has = q < R.p_vc.size() && (size_t) t < R.p_vc[q].size();
        f << '\t'; if (has) wr(f, R.p_vc[q][t]);  else f << "NA";
        f << '\t'; if (has) wr(f, R.p_acv[q][t]); else f << "NA";
        for (size_t k = 0; k < ctx.alphas.size(); ++k) {
            f << '\t';
            if (has && k < R.p_acatv[q][t].size()) wr(f, R.p_acatv[q][t][k]); else f << "NA";
        }
    };
    auto wr_alpha = [&](std::ofstream& f, const ChunkResultA& R, size_t q, int t) {
        if (!multi_a) return;
        for (size_t k = 0; k < ctx.alphas.size(); ++k) {
            f << '\t';
            if (q < R.p_alpha.size() && (size_t) t < R.p_alpha[q].size() && k < R.p_alpha[q][t].size())
                wr(f, R.p_alpha[q][t][k]);
            else f << "NA";
        }
    };
    long n1 = 0, n2 = 0;
    if (!pr.out_file.empty() && need_chunks) {
        std::ofstream f1(pr.out_file.c_str());
        if (!f1.is_open()) stop("Could not open out_file: " + pr.out_file);
        f1 << "chr\tstart\tend\twindow\tcategory\tm_cat\tm_flank\tphenotype\tvg\tse_vg\th2\tvg_flank\tvg_env";
        if (pr.spa)   f1 << "\tp_spa\tspa_used" << alpha_hdr << acv_hdr;
        if (pr.coher) f1 << "\tvg_t1\tvg_t2";
        if (pr.do_windows) {
            f1 << "\tvg_adj\th2_adj";
            if (pr.coher) f1 << "\tvg_t1_adj\tvg_t2_adj";
        }
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
                    if (pr.spa) { f1 << '\t'; wr(f1, R.p_spa[q][t]); f1 << '\t' << R.spa_used[q][t];
                                  wr_alpha(f1, R, q, t); wr_acv(f1, R, q, t); }
                    if (pr.coher) { f1 << '\t'; wr(f1, R.vg_t1[q][t]); f1 << '\t'; wr(f1, R.vg_t2[q][t]); }
                    if (pr.do_windows) {
                        const double a = adj[k][q][t];
                        f1 << '\t'; wr(f1, a); f1 << '\t'; wr(f1, (!std::isnan(a) && h2den(t) > 0) ? a / h2den(t) : NA_REAL);
                        if (pr.coher) {
                            f1 << '\t'; wr(f1, adj1[k][q][t]);
                            f1 << '\t'; wr(f1, adj2[k][q][t]);
                        }
                    }
                    f1 << '\n';
                }
        }
        f1.close();
        if (pr.do_windows) {
            std::ofstream f2((pr.out_file + ".windows").c_str());
            if (!f2.is_open()) stop("Could not open " + pr.out_file + ".windows");
            f2 << "chr\tstart\tend\twindow\tn_chunks\tcategory\tm_cat\tm_flank\tprs_variant\tphenotype\tvg\th2\tvg_flank\tvg_env";
            if (pr.coher) f2 << "\tvg_t1\tvg_t2";
            f2 << "\n";
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
                        f2 << '\t'; wr(f2, W.vg_flank[t]); f2 << '\t'; wr(f2, W.vg_env[t]);
                        if (pr.coher) { f2 << '\t'; wr(f2, W.vg_t1[kw][t]); f2 << '\t'; wr(f2, W.vg_t2[kw][t]); }
                        f2 << '\n';
                    }
                }
            f2.close();
        }
    }

    // ---- GENE WRITE ----------------------------------------------------------
    // Per gene: one row per (block, category), then one ACAT row combining every
    // p-value of that gene. vg_adj applies the window's per-unit-weight offset
    // delta_w to the gene's own cat_w -- the same rule the chunks get, so gene
    // and chunk heritabilities are calibrated by one definition rather than two.
    long ng = 0;
    if (!pr.out_file.empty() && has_genes) {
        std::ofstream f3((pr.out_file + ".genes").c_str());
        if (!f3.is_open()) stop("Could not open " + pr.out_file + ".genes");
        f3 << "gene\tchr\tstart\tend\tblock\tn_blocks\tcategory\tm_cat\tphenotype\tvg\tse_vg\th2\tvg_env";
        if (pr.spa)        f3 << "\tp_spa\tspa_used" << alpha_hdr << acv_hdr;
        if (pr.coher)      f3 << "\tvg_t1\tvg_t2";          // univariate estimates of both traits (for rg)
        if (pr.do_windows) f3 << "\tvg_adj\th2_adj";
        if (pr.cond_common) f3 << "\tn_cond\tcond_snps";
        f3 << "\n";
        // Window overlap per gene block, computed ONCE. Windows tile each
        // chromosome, so a block's SNP range intersects a short contiguous run
        // of them; doing this inside the trait loop would be O(windows) a row.
        std::vector< std::vector< std::pair<size_t,int> > > blk_win(gj.size());
        if (pr.do_windows)
            for (size_t b = 0; b < gj.size(); ++b) {
                if (!okg[b]) continue;
                for (size_t w = 0; w < wj.size(); ++w) {
                    if (wj[w].ci != gj[b].ci) continue;
                    const int lo = std::max(wj[w].a, gj[b].a), hi = std::min(wj[w].b, gj[b].b);
                    if (hi > lo) blk_win[b].push_back(std::make_pair(w, hi - lo));
                }
            }
        for (size_t gi = 0; gi < grec.size(); ++gi) {
            const GeneRec& R = grec[gi];
            const int nb = (int) (R.j1 - R.j0);
            bool any = false;
            for (size_t b = R.j0; b < R.j1; ++b) if (okg[b]) any = true;
            if (!any) continue;
            ++ng;
            for (int t = 0; t < NOUT; ++t) {
                std::vector<double> pv, pw;
                for (size_t b = R.j0; b < R.j1; ++b) {
                    if (!okg[b]) continue;
                    const ChunkResultA& B = rg[b];
                    for (size_t q = 0; q < B.cat_m.size(); ++q) {
                        // window-weighted correction for this block x category
                        double dnum = 0.0, dden = 0.0;
                        for (size_t z = 0; z < blk_win[b].size(); ++z) {
                            const double d = delta[blk_win[b][z].first][B.cat_id[q]][t];
                            if (std::isnan(d)) continue;
                            dnum += d * (double) blk_win[b][z].second;
                            dden += (double) blk_win[b][z].second;
                        }
                        const double a = (dden > 0.0)
                            ? B.vg[q][t] + (dnum / dden) * B.cat_w[q] : NA_REAL;
                        f3 << R.name << '\t' << R.chr << '\t' << R.start << '\t' << R.end << '\t'
                           << (int) (b - R.j0) + 1 << '\t' << nb << '\t'
                           << B.cat_name[q] << '\t' << B.cat_m[q] << '\t' << label(t) << '\t';
                        wr(f3, B.vg[q][t]); f3 << '\t'; wr(f3, B.se_vg[q][t]); f3 << '\t';
                        wr(f3, h2den(t) > 0 ? B.vg[q][t] / h2den(t) : NA_REAL);
                        f3 << '\t'; wr(f3, B.vg_env[t]);
                        if (pr.spa) {
                            f3 << '\t'; wr(f3, B.p_spa[q][t]); f3 << '\t' << B.spa_used[q][t];
                            wr_alpha(f3, B, q, t);
                            wr_acv(f3, B, q, t);
                            pv.push_back(B.p_spa[q][t]); pw.push_back((double) B.cat_m[q]);
                        }
                        if (pr.coher) { f3 << '\t'; wr(f3, B.vg_t1[q][t]); f3 << '\t'; wr(f3, B.vg_t2[q][t]); }
                        if (pr.do_windows) {
                            f3 << '\t'; wr(f3, a);
                            f3 << '\t'; wr(f3, (!std::isnan(a) && h2den(t) > 0) ? a / h2den(t) : NA_REAL);
                        }
                        if (pr.cond_common) f3 << '\t' << B.n_cond << '\t' << (B.cond_snps.empty() ? std::string("-") : B.cond_snps);
                        f3 << '\n';
                    }
                }
                if (!pr.spa || pv.empty()) continue;
                int mtot = 0; for (size_t z = 0; z < pw.size(); ++z) mtot += (int) pw[z];
                f3 << R.name << '\t' << R.chr << '\t' << R.start << '\t' << R.end << '\t'
                   << "NA\t" << nb << "\tACAT\t" << mtot << '\t' << label(t)
                   << "\tNA\tNA\tNA\tNA\t";
                wr(f3, acat_combine(pv, pw)); f3 << "\tNA";
                if (multi_a) for (size_t k = 0; k < ctx.alphas.size(); ++k) f3 << "\tNA";
                if (acv) for (size_t k = 0; k < ctx.alphas.size() + 2; ++k) f3 << "\tNA";
                if (pr.coher) f3 << "\tNA\tNA";
                if (pr.do_windows) f3 << "\tNA\tNA";
                if (pr.cond_common) f3 << "\tNA\tNA";
                f3 << '\n';
            }
        }
        f3.close();
    }

    if (need_chunks) {
        long skipped = 0; for (size_t k = 0; k < cj.size(); ++k) if (!ok1[k]) ++skipped;
        Rcout << "Chunks tested: " << n1 << (skipped ? "  (skipped " + std::to_string(skipped) + ")" : "") << "\n";
    }
    if (pr.do_windows) Rcout << "Window x category estimates: " << n2 << "\n";
    if (has_genes)     Rcout << "Genes tested: " << ng << " of " << grec.size() << "\n";
    return List::create(_["n_chunks"] = (double) n1, _["n_windows"] = (double) wj.size(),
                        _["n_genes"] = (double) ng,
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
                         int max_window_chunks = 8,
                         Rcpp::Nullable<Rcpp::CharacterMatrix> genes = R_NilValue,
                         int max_gene_snps = 512,
                         int min_gene_snps = 2,
                         Rcpp::NumericVector alpha = Rcpp::NumericVector::create(-1.0),
                         Rcpp::Nullable<Rcpp::CharacterVector> kvik_summaries = R_NilValue,
                         double acatv_mac = 10,
                         int collapse_mac = 0,
                         bool h2_without_prs = true,
                         bool one_sided = true,
                         double exact_ec = 10,
                         Rcpp::Nullable<Rcpp::NumericMatrix> covariates = R_NilValue,
                         double cov_df = NA_REAL,
                         bool SPA = true,
                         double spa_pval_threshold = NA_REAL,
                         bool binary = false,
                         Rcpp::Nullable<Rcpp::NumericMatrix> binary_raw = R_NilValue,
                         Rcpp::Nullable<Rcpp::NumericMatrix> binary_pi = R_NilValue,
                         bool coher = false,
                         bool cond_common = false,
                         double cond_p = 1e-4,
                         int cond_max = 10,
                         double cond_maf = 0.01,
                         SEXP chr = R_NilValue,
                         std::string out_file = "",
                         int batch_size = 64,
                         int n_threads = 0) {
    if (chunk_size < 1) stop("chunk_size must be >= 1");
    if (window_bp < chunk_size) stop("window_bp is smaller than one chunk");
    if (max_window_chunks < 1) stop("max_window_chunks must be >= 1");
    if (max_gene_snps < 1) stop("max_gene_snps must be >= 1");
    if (min_gene_snps < 1) stop("min_gene_snps must be >= 1");
    Rcout << "StratGWAS C++ build: " << __DATE__ << " " << __TIME__ << "   (p-values: one-sided VC / two-sided coher [conditional REML-GLS test; binary exact+lattice, exact Bernoulli rest, mode-centred lattice; quantitative: joint nuisance covariates + delta method for trait-1 components], cond_common (lead common variants as fixed effects; KVIK or correlation p), exact binary SPA + lattice, binary_pi, ACAT-V via KVIK, h2 without PRS, collapse_mac, no cumulant fallback (binary: NA if exact fails), coher vg_t1/vg_t2 in .genes)\n";
    // alpha may be a vector: alphas[0] is the primary alpha (vg, se, h2, the
    // window adjustment); with SPA and >1 value, p_spa is the ACAT over alphas.
    std::vector<double> alphas;
    for (int k = 0; k < (int) alpha.size(); ++k) {
        const double ak = alpha[k];
        if (!std::isfinite(ak)) stop("alpha must be finite");
        bool dup = false;
        for (size_t q = 0; q < alphas.size(); ++q) if (alphas[q] == ak) dup = true;
        if (!dup) alphas.push_back(ak);
    }
    if (alphas.empty()) stop("alpha must have at least one value");
    ChunkContext ctx = setup_context(filename, pheno_mat, alphas[0], covariates, annotation, annot_names);
    ctx.alphas = alphas;
    if (collapse_mac < 0) stop("collapse_mac must be >= 0");
    ctx.collapse_mac = collapse_mac;
    if (collapse_mac > 0)
        Rcout << "Collapsing: per tested unit and category, variants with MAC <= " << collapse_mac
              << " are summed into one pseudo-marker (ACAT-V still uses them individually).\n";
    if (alphas.size() > 1) {
        if (!SPA) Rcpp::warning("alpha has several values but SPA = FALSE: only alpha[1] is used");
        else {
            Rcout << "Multi-alpha: association p-values are the ACAT combination over alpha = {";
            for (size_t k = 0; k < alphas.size(); ++k) Rcout << (k ? ", " : "") << alphas[k];
            Rcout << "}; vg / se / h2 use alpha = " << alphas[0] << ".\n";
        }
    }
    const int P = ctx.n_pheno;

    // ---- gene definitions: name, chr, start, end (extra columns ignored) ------
    // Rows with a missing name or a non-numeric / missing coordinate are dropped
    // here rather than in the driver, so the count reported is what was usable.
    if (genes.isNotNull()) {
        Rcpp::CharacterMatrix gm(genes.get());
        if (gm.ncol() < 4) stop("genes must have at least 4 columns: name, chr, start, end");
        long bad = 0;
        // as.matrix() on a mixed data.frame pads via format(), so trim everything
        auto trim = [](const std::string& s) -> std::string {
            const size_t a = s.find_first_not_of(" \t\r\n");
            if (a == std::string::npos) return std::string();
            return s.substr(a, s.find_last_not_of(" \t\r\n") - a + 1); };
        for (int i = 0; i < gm.nrow(); ++i) {
            const std::string nm = trim(as<std::string>(gm(i, 0)));
            const std::string cs = as<std::string>(gm(i, 1));
            const std::string ss = as<std::string>(gm(i, 2));
            const std::string es = as<std::string>(gm(i, 3));
            if (nm.empty() || nm == "-" || nm == "NA" || cs.empty() || cs == "NA" ||
                ss.empty() || ss == "NA" || es.empty() || es == "NA") { ++bad; continue; }
            char* endp = 0;
            const long s = std::strtol(ss.c_str(), &endp, 10);
            if (endp == ss.c_str() || *endp != '\0') { ++bad; continue; }
            const long e = std::strtol(es.c_str(), &endp, 10);
            if (endp == es.c_str() || *endp != '\0') { ++bad; continue; }
            if (e < s) { ++bad; continue; }
            ctx.gene_name.push_back(nm); ctx.gene_chr.push_back(cs);
            ctx.gene_start.push_back(s); ctx.gene_end.push_back(e);
        }
        if (ctx.gene_name.empty()) stop("No usable rows in `genes`");
        Rcout << "Gene definitions: " << ctx.gene_name.size() << " read";
        if (bad) Rcout << ", " << bad << " dropped (missing name or coordinate)";
        Rcout << ". NOTE: these coordinates must be on the same genome build as the .bim.\n";
    }

    ChunkParams pr;
    pr.chunk_size = chunk_size; pr.window_bp = (long) window_bp; pr.do_windows = do_windows;
    pr.max_window_chunks = max_window_chunks;
    pr.max_gene_snps = max_gene_snps; pr.min_gene_snps = min_gene_snps;
    pr.spa = SPA; pr.binary = binary; pr.acatv_mac = acatv_mac; pr.h2_raw = false;
    pr.one_sided = one_sided;
    if (SPA) Rcout << (one_sided ? "Variance-component p-values: one-sided (upper tail)."
                                 : "Variance-component p-values: two-sided (from the median).")
                   << (coher ? " Co-heritability p-values are always two-sided." : "") << "\n";
    if (!ISNAN(spa_pval_threshold))
        Rcout << "NOTE: spa_pval_threshold is no longer used -- every test now uses the saddlepoint.\n";
    pr.cov_df = ISNAN(cov_df) ? ((int) ctx.covZ.cols() + 1) : (int) cov_df;
    if (pr.cov_df < 0) pr.cov_df = 0;
    if (pr.cov_df >= ctx.n_inds) stop("cov_df must be smaller than the sample size");
    Rcout << "Environment moment: T(env,env) = n - " << pr.cov_df << "\n";
    pr.out_file = out_file; pr.batch_size = batch_size; pr.n_threads = n_threads;
    pr.cond_common = cond_common; pr.cond_p = cond_p; pr.cond_max = cond_max; pr.cond_maf = cond_maf;
    if (cond_common) {
        if (kvik_summaries.isNull())
            Rcout << "cond_common: lead common variants selected by their correlation with the phenotypes "
                  << "(z = r sqrt(n); no kvik_summaries given), p < " << cond_p << " for any phenotype, at most "
                  << cond_max << " per gene block, MAF >= " << cond_maf << ".\n";
        if (!(cond_p > 0.0 && cond_p < 1.0)) stop("cond_p must be in (0, 1)");
        if (cond_max < 1) stop("cond_max must be >= 1");
    }
    pr.use_acatv = SPA && !coher && !cond_common && kvik_summaries.isNotNull();

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
    pr.h2_raw = (G > 0) && h2_without_prs;
    if (pr.h2_raw)
        Rcout << "LOCO PRS used for p-values only: vg / se_vg / h2 of chunks and genes are estimated on y itself\n"
              << "  (an in-sample PRS shrinks heritability estimates). Use prs_mask = 0 to estimate windows on y as well.\n";
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
    // ---- raw 0/1 phenotypes: binary co-heritability (required) and the exact ------
    // ---- per-variant saddlepoint for single binary traits (optional) ----------------
    if (binary && binary_raw.isNull())
        stop("binary = TRUE needs binary_raw (the original 0/1 phenotypes, one column per trait); sparrow() passes it automatically");
    if (binary && binary_raw.isNotNull()) {
        Rcpp::NumericMatrix braw(binary_raw.get());
        if (braw.ncol() != P) stop("binary_raw must have one column per phenotype");
        if (braw.nrow() != ctx.Y.rows() && braw.nrow() < 1) stop("binary_raw has no rows");
        const int n = ctx.n_inds;
        ctx.braw.assign(P, std::vector<unsigned char>(n, 0));
        ctx.prev.assign(P, 0.0); ctx.bdelta.assign(P, 0.0);
        ctx.bpi.assign(P, std::vector<float>(n, 0.f)); ctx.ncase.assign(P, 0.0);
        for (int t = 0; t < P; ++t) {
            // y = delta * (b - mu): mu = covariate fit of the 0/1 trait. delta from y'y = delta b'y
            // (exact for a least-squares residual of b, scaled or not).
            double ym = 0.0;
            for (int i = 0; i < n; ++i) ym += ctx.Y(i, t);
            ym /= n;
            double sumb = 0.0, bty = 0.0, yy = 0.0;
            for (int i = 0; i < n; ++i) {
                const double v = braw(ctx.pheno_keep[i], t);
                if (!(v == 0.0 || v == 1.0)) stop("binary_raw must contain only 0 and 1 (column %d)", t + 1);
                const double yc = ctx.Y(i, t) - ym;
                ctx.braw[t][i] = (unsigned char) v; sumb += v; bty += v * yc; yy += yc * yc;
            }
            ctx.prev[t] = sumb / n; ctx.ncase[t] = sumb;
            if (!(bty > 0.0)) stop("binary_raw column %d does not match the phenotype", t + 1);
            ctx.bdelta[t] = yy / bty;
            // null P(case): the logistic null model when given (binary_pi, recommended: the
            // linear fit below is clipped at 0 and makes the exact saddlepoint conservative
            // when covariates predict the trait); otherwise the linear covariate fit of b.
            for (int i = 0; i < n; ++i) {
                double mu = (double) ctx.braw[t][i] - (ctx.Y(i, t) - ym) / ctx.bdelta[t];
                if (mu < 1e-8) mu = 1e-8; if (mu > 1.0 - 1e-8) mu = 1.0 - 1e-8;
                ctx.bpi[t][i] = (float) mu;
            }
        }
        if (binary_pi.isNotNull()) {
            Rcpp::NumericMatrix bp(binary_pi.get());
            if (bp.ncol() != P || bp.nrow() != braw.nrow()) stop("binary_pi must have the same dimensions as binary_raw");
            for (int t = 0; t < P; ++t)
                for (int i = 0; i < n; ++i) {
                    double mu = bp(ctx.pheno_keep[i], t);
                    if (!std::isfinite(mu) || mu < 0.0 || mu > 1.0) stop("binary_pi must contain probabilities (column %d)", t + 1);
                    ctx.bpi[t][i] = (float) std::min(std::max(mu, 1e-8), 1.0 - 1e-8);
                }
        }
        if (coher) {
            Rcout << "Binary co-heritability: exact conditional Bernoulli saddlepoint; prevalences:";
            for (int t = 0; t < P; ++t) Rcout << " " << ctx.prev[t];
            Rcout << "\n";
        }
        if (!coher && SPA) {
            ctx.bin_exact = true; ctx.car_max = 5000; ctx.ec_max = exact_ec;
            Rcout << "Binary phenotype: exact per-variant saddlepoint (spa_used = 5): exact case-count "
                  << "distribution for variants with <= " << exact_ec << " expected case carriers, "
                  << "Gaussian for the rest, lattice convolution when p < 0.1; null P(case) from "
                  << (binary_pi.isNotNull() ? "the logistic null model" : "the linear covariate fit (pass binary_pi)") << ".\n";
        }
    }

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
    // ---- LDAK-KVIK single-variant summaries -> ACAT-V --------------------------
    if (kvik_summaries.isNotNull()) {
        Rcpp::CharacterVector kv(kvik_summaries.get());
        if (!pr.use_acatv && !pr.cond_common) {
            if (coher) Rcpp::warning("kvik_summaries is ignored with coher = TRUE (ACAT-V is univariate)");
            else Rcpp::warning("kvik_summaries is ignored with SPA = FALSE");
        } else {
            if (kv.size() != P) stop("kvik_summaries must give one file per phenotype (column of pheno)");
            std::unordered_map<std::string, int> idx;
            idx.reserve((size_t) ctx.n_snps * 2);
            long dup = 0;
            for (int j = 0; j < ctx.n_snps; ++j)
                if (!idx.insert(std::make_pair(ctx.bim.id[j], j)).second) ++dup;
            if (dup) Rcout << "NOTE: " << dup << " duplicated SNP names in the .bim; KVIK results are matched to the first.\n";
            ctx.kz.resize(P); ctx.kmac.resize(P);
            for (int t = 0; t < P; ++t) kvik_load(as<std::string>(kv[t]), ctx, t, idx);
            if (pr.cond_common)
                Rcout << "cond_common: each gene block is tested conditional on its lead common variants "
                      << "(MAF >= " << pr.cond_maf << ", KVIK p < " << pr.cond_p << " for any trait, at most "
                      << pr.cond_max << ", r^2 < 0.9 with those already selected), projected out of genotypes "
                      << "and phenotypes" << (pr.binary ? " (binary: null P(case) updated)" : "") << ".\n"
                      << "  ACAT-V (marginal single-variant tests) is not used with cond_common.\n";
            if (pr.use_acatv) Rcout << "ACAT-V: KVIK single-variant tests, weights [2f(1-f)]^(1+alpha) for each alpha, "
                  << "variants with MAC < " << acatv_mac << " collapsed into a burden.\n"
                  << "  p_spa = ACAT over " << 2 * ctx.alphas.size() << " p-values ("
                  << ctx.alphas.size() << " variance-component + " << ctx.alphas.size() << " ACAT-V).\n";
        }
    }
    return run_all(ctx, pr);
}
