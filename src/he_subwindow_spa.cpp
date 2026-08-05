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
// CHUNK-BASED HE + saddlepoint association testing.
//
// The unit of analysis is a fixed-size CHUNK of consecutive SNPs (default
// 256), not a bp window. That single change is what makes the method
// tractable, because it puts a HARD BOUND on the only quantity that matters
// for cost.
//
// WHY THE bp-WINDOW DESIGNS FAILED. The saddlepoint calculation needs, per
// tested component, a Cholesky (O(K^3)) and a symmetric eigensolve (~10 K^3)
// on a K x K matrix, where K is the TOTAL column count over all components.
// With bp-defined windows K is set by the DATA, not by the user: a dense 1 Mb
// exome cell holds ~8547 SNPs, so
//     K = 8547 + 406 (common) = 8953
//     -> 641 MB per K x K matrix, ~12 of them live  =>  ~7.7 GB per thread
//     -> ~885 s per tested sub-window  =>  ~49 h per cell.
// No amount of algebraic reuse rescues that; the exponent is the problem.
//
// WITH CHUNKS, K IS CHOSEN BY THE USER:
//     K = chunk_size * (1 + 2 * flank_chunks) + m_common
// e.g. 256 * 3 + ~20 = ~788, giving K^3 ~ 4.9e8 -- roughly 1500x cheaper per
// test than K = 8953, and ~5 MB per matrix instead of 641 MB. Cost becomes
// predictable and independent of local SNP density.
//
// MODEL for the chunk at SNP indices [a, b):
//     { K_chunk (target), K_flank (the flank_chunks blocks either side,
//       combined), [K_common], sigma_e I }
// The flanking chunks play the role the 1 Mb flanks played -- absorbing local
// LD and background -- but at bounded size.
//
// TWO FIXES CARRIED OVER FROM THE EARLIER FILES:
//   * The Wald pre-screen is applied BEFORE the eigensolve, not after. Since
//     the explicit eigenvalues are those of the symmetric matrix Asym,
//         Var(sigma_hat) = 2 * sum lambda^2 = 2 * ||Asym||_F^2 ,
//     which costs O(K^2) given Asym -- no eigensolve needed. In the previous
//     file the screen sat AFTER the ~10 K^3 eigensolve and therefore saved
//     nothing at all. This ordering is what makes it a real ~10x saving.
//   * G is accumulated in float from near-collinear ultra-rare columns, so the
//     Cholesky needs a real ridge (1e-4 relative, escalating), not 1e-10.
//
// Everything is in an anonymous namespace (internal linkage) so it cannot
// collide with the identically-named helpers in the other .cpp files.
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

// Read a cell by bp range [cell_start, cell_start + W) within chromosome index
// range [lo, hi] of `bim`. Used for the COMMON fileset (matched to a WES cell's
// bp span). Thread-safe; false on I/O failure.
static bool read_cell(const BedReader& br, const BimInfo& bim,
                      const std::vector<int>& keep,
                      int lo, int hi, long cell_start, long W, Cell& cell) {
    cell = Cell(); cell.i0 = -1;
    if (lo < 0 || hi < lo || W <= 0) { cell.X = GenoMat(keep.size(), 0); return true; }
    int i0 = lower_index(bim.bp, lo, hi, cell_start);
    int i1 = lower_index(bim.bp, lo, hi, cell_start + W);
    int m = i1 - i0;
    if (m <= 0) { cell.X = GenoMat(keep.size(), 0); return true; }
    GenoMat X;
    if (!br.read_raw(keep, i0, i1, X)) return false;
    standardize_capture_maf(X, cell.maf);
    cell.X = std::move(X);
    cell.i0 = i0;
    return true;
}

// Partition a chromosome's WES SNP indices [c_lo, c_hi] into consecutive CELLS.
// A cell grows from its start until BOTH (span >= W bp) AND (>= min_snps SNPs),
// or the chromosome end is reached. A trailing remainder too small to form a
// full cell on its own is absorbed into the last cell (so no tiny tail window).
// Returns the [first,last] WES index of each cell.
static std::vector< std::pair<int,int> > build_cells(const std::vector<long>& bp,
        int c_lo, int c_hi, long W, int min_snps) {
    std::vector< std::pair<int,int> > cells;
    if (c_hi < c_lo) return cells;
    int s = c_lo;
    while (s <= c_hi) {
        int e = s;
        while (e < c_hi) {
            long span = bp[e] - bp[s] + 1;          // inclusive bp span so far
            int  cnt  = e - s + 1;
            if (span >= W && cnt >= min_snps) break;
            ++e;
        }
        // Absorb a trailing remainder that could not itself satisfy the rules.
        int  rem_cnt  = c_hi - e;                    // SNPs after e
        long rem_span = (e < c_hi) ? (bp[c_hi] - bp[e + 1] + 1) : 0;
        if (rem_cnt > 0 && (rem_cnt < min_snps || rem_span < W)) e = c_hi;
        cells.push_back(std::make_pair(s, e));
        s = e + 1;
    }
    return cells;
}


// ===========================================================================
// Context
// ===========================================================================
struct ChunkContext {
    std::string wes_prefix; int wes_n_total, wes_n_snps; BimInfo wes_bim;
    std::vector<int> geno_keep;
    Eigen::MatrixXd Y;
    CharacterVector trait_names;
    int n_inds, n_pheno;
    std::vector<std::string> chr_order;
    std::vector<int> chr_lo, chr_hi;
    bool use_common; std::string common_prefix;
    int common_n_total, common_n_snps; BimInfo common_bim;
    std::vector<int> common_keep;
    double alpha, alpha_common;
    std::vector<float> snp_weights;
    GenoMat covZ, covM;
    // Optional annotation: EVERY column is a tested functional category.
    // snp_cat[j] = 0..n_annot_cat-1 (a tested category) or -1 (flank/background).
    bool use_annot; int n_annot_cat;
    std::vector<std::string> annot_names;
    std::vector<int> snp_cat;
    // Thread-safe genotype readers (replace readBedBlock in the hot path).
    BedReader wes, com;
    // Precomputed per-chromosome index range in the COMMON .bim. The old code
    // recomputed this with a full linear scan over every common SNP for every
    // chunk (~1.2e10 string compares at 1M common SNPs x 11.7k chunks), inside
    // the serial section.
    std::map<std::string, std::pair<int,int> > common_range;
};

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

static ChunkContext setup_chunk_context(const std::string& filename, const SEXP pheno_mat,
                                    double alpha, double alpha_common,
                                    Rcpp::Nullable<Rcpp::String> common_filename,
                                    Rcpp::Nullable<Rcpp::NumericVector> weights,
                                    Rcpp::Nullable<Rcpp::NumericMatrix> covariates,
                                    Rcpp::Nullable<Rcpp::IntegerMatrix> annotation,
                                    Rcpp::Nullable<Rcpp::CharacterVector> annot_names) {
    ChunkContext ctx;
    ctx.use_annot = false; ctx.n_annot_cat = 0;
    ctx.wes_prefix = filename;
    ctx.wes_n_snps = count_lines(filename + ".bim");
    List fam = read_fam_file(filename);
    CharacterVector geno_iid = fam["iid"];
    ctx.wes_n_total = geno_iid.size();
    ctx.alpha = alpha; ctx.alpha_common = alpha_common;
    Rcout << "WES individuals (.fam): " << ctx.wes_n_total << "\n";
    Rcout << "WES SNPs (.bim): "        << ctx.wes_n_snps  << "\n";
    ctx.wes_bim = read_bim_positions(filename + ".bim");
    if (ctx.wes_bim.n_snps != ctx.wes_n_snps) stop("Mismatch between .bim line count and parsed positions");
    calibrate_bed_reader(ctx.wes, filename, ctx.wes_n_total, ctx.wes_n_snps, "WES");

    if (weights.isNotNull()) {
        Rcpp::NumericVector wv(weights.get());
        if (wv.size() != ctx.wes_n_snps) stop("weights must have one entry per SNP in the WES .bim file");
        ctx.snp_weights.resize(ctx.wes_n_snps);
        for (int j = 0; j < ctx.wes_n_snps; ++j) {
            double w = wv[j]; if (ISNAN(w) || w < 0.0) w = 0.0;
            ctx.snp_weights[j] = (float) w;
        }
        Rcout << "Using per-SNP LD weights\n";
    }

    // Annotation: EVERY column is a tested functional category -- there is no
    // flank-flag column any more. The background/flanking component is built
    // purely from the ADJACENT CHUNKS (flank_chunks) plus, if supplied, the
    // common SNPs. A SNP whose annotation row is all zero belongs to no category
    // and simply joins that background, so no genotype is silently dropped.
    if (annotation.isNotNull()) {
        Rcpp::IntegerMatrix am(annotation.get());
        if (am.nrow() != ctx.wes_n_snps)
            stop("annotation must have one row per SNP in the WES .bim file");
        int ncol = am.ncol();
        if (ncol < 1) stop("annotation needs at least one column; every column is a tested category");
        ctx.n_annot_cat = ncol;
        ctx.snp_cat.assign(ctx.wes_n_snps, -1);
        for (int j = 0; j < ctx.wes_n_snps; ++j)
            for (int c = 0; c < ncol; ++c) if (am(j, c) == 1) { ctx.snp_cat[j] = c; break; }
        if (annot_names.isNotNull()) {
            Rcpp::CharacterVector nm(annot_names.get());
            if (nm.size() != ctx.n_annot_cat)
                stop("annot_names length must equal the number of annotation columns "
                     "(every column is a category)");
            for (int c = 0; c < ctx.n_annot_cat; ++c) ctx.annot_names.push_back(as<std::string>(nm[c]));
        } else {
            for (int c = 0; c < ctx.n_annot_cat; ++c) ctx.annot_names.push_back("cat" + std::to_string(c + 1));
        }
        ctx.use_annot = true;
        Rcout << "Annotation: " << ctx.n_annot_cat
              << " categories (every column); an SPA p-value is reported per category.\n"
              << "  Background = adjacent chunks + common SNPs (not taken from the annotation).\n";
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
    std::vector<int> pheno_keep;
    for (int i = 0; i < match_idx.size(); ++i)
        if (match_idx[i] != NA_INTEGER) { ctx.geno_keep.push_back(i); pheno_keep.push_back(match_idx[i] - 1); }
    if (ctx.geno_keep.empty()) stop("No overlapping individuals between genotype and phenotype");
    ctx.n_inds = (int) ctx.geno_keep.size();
    Rcout << "Individuals with complete data: " << ctx.n_inds << "\n";

    ctx.Y.resize(ctx.n_inds, ctx.n_pheno);
    for (int i = 0; i < ctx.n_inds; ++i)
        for (int j = 0; j < ctx.n_pheno; ++j) ctx.Y(i, j) = pheno(pheno_keep[i], j);

    if (covariates.isNotNull()) {
        Rcpp::NumericMatrix cv(covariates.get());
        if (cv.nrow() != pheno.nrow())
            stop("covariates must have one row per row of the phenotype matrix");
        int qc = cv.ncol();
        if (qc > 0) {
            Eigen::MatrixXd Zd(ctx.n_inds, qc);
            for (int i = 0; i < ctx.n_inds; ++i)
                for (int j = 0; j < qc; ++j) { double v = cv(pheno_keep[i], j); Zd(i, j) = ISNAN(v) ? 0.0 : v; }
            for (int j = 0; j < qc; ++j) Zd.col(j).array() -= Zd.col(j).mean();
            Eigen::MatrixXd ZtZ = Zd.transpose() * Zd;
            Eigen::MatrixXd Md = ZtZ.completeOrthogonalDecomposition().pseudoInverse() * Zd.transpose();
            ctx.covZ = Zd.cast<float>(); ctx.covM = Md.cast<float>();
            Rcout << "Regressing " << qc << " covariates out of genotypes\n";
        }
    }

    CharacterVector analysis_iid(ctx.n_inds);
    for (int i = 0; i < ctx.n_inds; ++i) analysis_iid[i] = geno_iid[ctx.geno_keep[i]];
    for (int j = 0; j < ctx.wes_n_snps; ++j) {
        if (ctx.chr_order.empty() || ctx.wes_bim.chr[j] != ctx.chr_order.back()) {
            ctx.chr_order.push_back(ctx.wes_bim.chr[j]);
            ctx.chr_lo.push_back(j); ctx.chr_hi.push_back(j);
        } else ctx.chr_hi.back() = j;
    }

    ctx.use_common = common_filename.isNotNull();
    ctx.common_n_total = ctx.common_n_snps = 0;
    if (ctx.use_common) {
        ctx.common_prefix = as<std::string>(common_filename.get());
        ctx.common_n_snps = count_lines(ctx.common_prefix + ".bim");
        List cfam = read_fam_file(ctx.common_prefix);
        CharacterVector common_iid = cfam["iid"];
        ctx.common_n_total = common_iid.size();
        ctx.common_bim = read_bim_positions(ctx.common_prefix + ".bim");
        Rcout << "Common SNPs (.bim): " << ctx.common_n_snps << "\n";
        calibrate_bed_reader(ctx.com, ctx.common_prefix, ctx.common_n_total, ctx.common_n_snps, "common");
        for (int j = 0; j < ctx.common_n_snps; ++j) {
            const std::string& c = ctx.common_bim.chr[j];
            std::map<std::string, std::pair<int,int> >::iterator it = ctx.common_range.find(c);
            if (it == ctx.common_range.end()) ctx.common_range[c] = std::make_pair(j, j);
            else it->second.second = j;
        }
        IntegerVector cmatch = match(analysis_iid, common_iid);
        ctx.common_keep.resize(ctx.n_inds);
        for (int i = 0; i < ctx.n_inds; ++i) {
            if (cmatch[i] == NA_INTEGER) stop("An analysis individual is missing from the common-SNP .fam file");
            ctx.common_keep[i] = cmatch[i] - 1;
        }
    }
    return ctx;
}

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
struct QuadSpaResult { double p; bool converged; };

static void quad_cgf(double t, const std::vector<double>& eig_explicit,
                     double eig_rep, long n_rep, double& K, double& K1, double& K2) {
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
        double nr = (double) n_rep;
        K  += -0.5 * nr * std::log(d);
        K1 += nr * eig_rep / d;
        K2 += nr * 2.0 * eig_rep * eig_rep / (d * d);
    }
}

static QuadSpaResult quad_spa_solve(
    double s_obs_in, const std::vector<double>& eig_in, double eig_rep_in, long n_rep,
    int max_iter = 100, double tol = 1e-8
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
        mean_raw += (double) n_rep * eig_rep_in;
        var_raw  += (double) n_rep * 2.0 * eig_rep_in * eig_rep_in;
    }
    // Only a genuinely degenerate (zero / non-finite) variance is unusable.
    if (!(var_raw > 0.0) || !std::isfinite(var_raw)) return res;

    const double sd_raw = std::sqrt(var_raw);
    const double scale  = 1.0 / sd_raw;

    std::vector<double> eig_explicit(eig_in.size());
    for (size_t j = 0; j < eig_in.size(); ++j) eig_explicit[j] = eig_in[j] * scale;
    const double eig_rep = eig_rep_in * scale;
    const double s_obs   = s_obs_in   * scale;

    double lam_min = eig_rep, lam_max = eig_rep;
    for (size_t j = 0; j < eig_explicit.size(); ++j) {
        double lam = eig_explicit[j];
        if (lam < lam_min) lam_min = lam;
        if (lam > lam_max) lam_max = lam;
    }
    // Dimensionless now: |lambda| <= 1/sqrt(2), so 1e-14 simply means
    // "no positive (negative) eigenvalue, hence no bound on t in that direction".
    const double lam_eps = 1e-14;
    double t_hi = (lam_max >  lam_eps) ? (1.0 / (2.0 * lam_max)) :  1e12;
    double t_lo = (lam_min < -lam_eps) ? (1.0 / (2.0 * lam_min)) : -1e12;
    double margin = 1e-6;
    double t_hi_safe = (lam_max >  lam_eps) ? t_hi * (1.0 - margin) : t_hi;
    double t_lo_safe = (lam_min < -lam_eps) ? t_lo * (1.0 - margin) : t_lo;

    double t = 0.0, K, K1, K2;
    quad_cgf(0.0, eig_explicit, eig_rep, n_rep, K, K1, K2);
    double mean0 = K1, var0 = K2;          // var0 == 1 up to rounding
    if (!(var0 > 0.0)) return res;

    bool converged = false;
    for (int it = 0; it < max_iter; ++it) {
        quad_cgf(t, eig_explicit, eig_rep, n_rep, K, K1, K2);
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

    quad_cgf(t, eig_explicit, eig_rep, n_rep, K, K1, K2);
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

    double p_one = (w >= 0)
        ? (1.0 - Phi_w) + phi_w * (1.0 / u - 1.0 / w)
        : Phi_w - phi_w * (1.0 / u - 1.0 / w);
    if (p_one < 0.0) p_one = 0.0;
    if (p_one > 1.0) p_one = 1.0;

    res.p = std::min(1.0, 2.0 * p_one);
    res.converged = true;
    return res;
}



// ===========================================================================
// One chunk's assembled region
// ===========================================================================
struct ChunkData {
    std::string chr; long start, end;
    int i0_target, m_t, m_f, m_c, K;   // target / flank / common column counts
    GenoMat V;                          // n x K = [target | flank | common]
    std::vector<double> cj;             // ||v_j||^2 per column
    // When flank_chunks = 0 the common SNPs ARE the background component, so they
    // occupy the flank block of V. Only affects how the columns are reported.
    bool flank_is_common;
    ChunkData() : i0_target(0), m_t(0), m_f(0), m_c(0), K(0), flank_is_common(false) {}
};

struct ChunkResult {
    std::string chr; long start, end;
    int m_t, m_f, m_c;
    std::vector<double> vg, se_vg, h2, p_spa;
    std::vector<int> spa_used;
};

// Components: 0 = target chunk, 1 = flank, 2 = common (if any), env = residual.
// Kernels use tr(K_a) = n, i.e. K_a = g_a S_A with g_a = n / tr(S_A), so
// T(a,env) = n exactly for every a.
static void test_chunk(const ChunkData& cd, const Eigen::MatrixXd& Y, const GenoMat& Yf,
                       const std::vector<double>& Vp, const std::vector<double>& yty,
                       bool spa, double spa_thresh,
                       ChunkResult& cr) {
    const int n = (int) Y.rows(), P = (int) Y.cols();
    const int K = cd.K, m_t = cd.m_t, m_f = cd.m_f, m_c = cd.m_c;
    const bool has_c = (m_c > 0);
    // Components present, in V-column order: target (always), flank (if any),
    // common (if any). Built dynamically so the model degrades gracefully when
    // there is no background at all -- flank_chunks = 0 with no common fileset
    // leaves { target, sigma_e I }, which is unconditioned but still estimable.
    std::vector<int> off;  off.push_back(0);  off.push_back(m_t);
    if (m_f > 0) off.push_back(m_t + m_f);
    if (m_c > 0) off.push_back(m_t + m_f + m_c);
    const int C = (int) off.size() - 1, env = C;

    cr.chr = cd.chr; cr.start = cd.start; cr.end = cd.end;
    // Report honestly: when flank_chunks = 0 the flank block holds the common SNPs.
    cr.m_t = m_t;
    cr.m_f = cd.flank_is_common ? 0   : m_f;
    cr.m_c = cd.flank_is_common ? m_f : m_c;
    cr.vg.assign(P, NA_REAL); cr.se_vg.assign(P, NA_REAL);
    cr.h2.assign(P, NA_REAL); cr.p_spa.assign(P, NA_REAL);
    cr.spa_used.assign(P, 0);
    if (m_t <= 0) return;

    // ---- Gram matrix (K x K, K is BOUNDED by chunk_size and flank_chunks) ---
    GenoMat Gf = GenoMat::Zero(K, K);
    Gf.selfadjointView<Eigen::Upper>().rankUpdate(cd.V.transpose());
    Eigen::MatrixXd G = GenoMat(Gf.selfadjointView<Eigen::Upper>()).cast<double>();
    Eigen::MatrixXd G2 = G.array().square();

    std::vector<double> g(C);
    for (int a = 0; a < C; ++a) {
        double tr = 0.0;
        for (int j = off[a]; j < off[a + 1]; ++j) tr += cd.cj[j];
        if (!(tr > 0.0)) return;
        g[a] = (double) n / tr;
    }

    // ---- moment matrix T:  tr(S_a S_b) = ||G[A,B]||_F^2 ---------------------
    Eigen::MatrixXd T = Eigen::MatrixXd::Zero(C + 1, C + 1);
    for (int a = 0; a < C; ++a)
        for (int b = a; b < C; ++b) {
            double v = g[a] * g[b] *
                G2.block(off[a], off[b], off[a+1]-off[a], off[b+1]-off[b]).sum();
            T(a,b) = v; T(b,a) = v;
        }
    for (int a = 0; a < C; ++a) { T(a,env) = (double) n; T(env,a) = (double) n; }
    T(env,env) = (double) n;
    Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> Tcod(T);
    Eigen::MatrixXd Tinv = Tcod.pseudoInverse();

    // ---- Cholesky of G, once per chunk -------------------------------------
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

    // ---- BLOCK OUTER PRODUCTS, once per chunk (see build_SA below) ---------
    // D0 and DM are block-constant -- they take only C distinct values, one per
    // component -- so  L^T D L = sum_a d_a * (L_a^T L_a)  with L_a the rows of L
    // belonging to component a. The W_a depend only on the chunk, not on the
    // trait, so hoisting them out of the trait loop turns two O(K^3) triple
    // products PER TRAIT into O(C K^2) linear combinations. Benchmarked 3.86x
    // at K=1168, C=3, P=20 (see the companion note).
    std::vector<Eigen::MatrixXd> W;
    bool use_wcache = spa && have_L && ((double) C * K * K * 8.0 <= 6e8);
    if (use_wcache) {
        W.resize(C);
        for (int a = 0; a < C; ++a) {
            const int ma = off[a + 1] - off[a];
            W[a].noalias() = L.middleRows(off[a], ma).transpose() * L.middleRows(off[a], ma);
        }
    }
    // Build S = sigma_env0 I + L^T D0 L  and  A = c_env I + L^T DM L.
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

    // ---- one GEMM for all traits, instead of a GEMV per trait --------------
    // V is n x K (467 MB at n=1e5, K=1168); the old per-trait
    // V.transpose()*Y.col(t) streamed all of it once PER TRAIT.
    Eigen::MatrixXd U = (cd.V.transpose() * Yf).cast<double>();      // K x P

    for (int t = 0; t < P; ++t) {
        const Eigen::VectorXd u = U.col(t);
        Eigen::VectorXd q(C + 1);
        for (int a = 0; a < C; ++a) {
            double qa = 0.0;
            for (int j = off[a]; j < off[a + 1]; ++j) qa += u[j] * u[j];
            q[a] = g[a] * qa;
        }
        q[env] = yty[t];                    // constant across chunks; hoisted

        Eigen::VectorXd sigma = Tcod.solve(q);
        cr.vg[t] = sigma[0];
        cr.h2[t] = (Vp[t] > 0) ? sigma[0] / Vp[t] : NA_REAL;
        if (!spa || !have_L) continue;

        // ---- plug-in null covariance ---------------------------------------
        Eigen::VectorXd s0v = sigma.head(C);
        for (int a = 0; a < C; ++a) if (s0v[a] < 0.0) s0v[a] = 0.0;
        s0v[0] = 0.0;
        double sigma_env0 = std::max(sigma[env], 1e-8 * (Vp[t] > 0 ? Vp[t] : 1.0));

        const double c_env = Tinv(0,env);
        for (int a = 0; a < C; ++a) { d0v[a] = s0v[a] * g[a]; dmv[a] = Tinv(0,a) * g[a]; }
        d0v[0] = 0.0;                                    // null the tested component
        build_SA();
        Smat.diagonal().array() += sigma_env0;           // S
        Amat.diagonal().array() += c_env;                // A = Q + c_env I

        // ---- WALD SCREEN, WITHOUT BUILDING Asym -----------------------------
        // The explicit eigenvalues are those of Asym = Cf^T A Cf with S = Cf Cf^T,
        // so with N = A S,
        //     Var(sigma_hat) = 2 ||Asym||_F^2 = 2 tr((A S)^2) = 2 sum_ij N_ij N_ji ,
        // which needs ONE K^3 product and an O(K^2) contraction -- no Cholesky of
        // S and no Cf products. The old code paid ~4.3 K^3 to build Asym before it
        // could screen at all, so the screen only ever saved the eigensolve
        // (~30%); now the expensive path is genuinely skipped for the ~90% of
        // tests that are nowhere near significance.
        double eig_rep = c_env * sigma_env0;
        long n_rep = (long) n - K;
        Eigen::MatrixXd N; N.noalias() = Amat * Smat;
        double var0 = 2.0 * (N.cwiseProduct(N.transpose())).sum();
        if (n_rep > 0) var0 += (double) n_rep * 2.0 * eig_rep * eig_rep;
        if (!(var0 > 0.0)) continue;
        cr.se_vg[t] = std::sqrt(var0);
        double p_wald = std::erfc(std::abs(sigma[0] / cr.se_vg[t]) / std::sqrt(2.0));
        if (p_wald >= spa_thresh) { cr.p_spa[t] = p_wald; cr.spa_used[t] = 0; continue; }

        // ---- only now build Asym and solve the eigenproblem ------------------
        Eigen::LLT<Eigen::MatrixXd> lltS(Smat);
        if (lltS.info() != Eigen::Success) { cr.p_spa[t] = p_wald; continue; }
        Eigen::MatrixXd Cf = lltS.matrixL();
        Eigen::MatrixXd Asym = Cf.transpose() * Amat * Cf;
        Asym = (0.5*(Asym + Asym.transpose())).eval();
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ses(Asym, Eigen::EigenvaluesOnly);
        if (ses.info() != Eigen::Success) { cr.p_spa[t] = p_wald; continue; }
        const Eigen::VectorXd& ev = ses.eigenvalues();
        std::vector<double> eig(ev.data(), ev.data() + K);
        QuadSpaResult qr = quad_spa_solve(sigma[0], eig, eig_rep, n_rep);
        if (qr.converged) { cr.p_spa[t] = qr.p; cr.spa_used[t] = 1; }
        else              { cr.p_spa[t] = p_wald; cr.spa_used[t] = 0; }
    }
}

// ===========================================================================
// Driver
// ===========================================================================
struct ChunkParams {
    int chunk_size, flank_chunks, min_chunk_snps;
    long common_bp; bool common_window_given; int max_common_snps;
    bool spa; double spa_thresh;
    std::string out_file; int batch_size, n_threads;
};

static void common_chr_range(const ChunkContext& ctx, const std::string& chr, int& lo, int& hi) {
    lo = hi = -1;
    if (!ctx.use_common) return;
    std::map<std::string, std::pair<int,int> >::const_iterator it = ctx.common_range.find(chr);
    if (it != ctx.common_range.end()) { lo = it->second.first; hi = it->second.second; }
}

// Find the [i0, i1) index range (within chromosome index range [lo, hi] of a
// BimInfo) of the `cap` SNPs whose bp is CLOSEST to `center`. Since bp is
// sorted, the cap nearest points always form a contiguous run, so this is a
// two-pointer expansion out from the insertion point -- O(log n) to locate it
// plus O(cap) to walk outward. Used when common_window is left unspecified:
// no bp window, no thinning, just the cap nearest common SNPs directly.
static std::pair<int,int> nearest_snp_range(const std::vector<long>& bp, int lo, int hi,
                                            long center, int cap) {
    if (lo < 0 || hi < lo || cap <= 0) return std::make_pair(0, 0);
    int ins = lower_index(bp, lo, hi, center);
    int L = ins - 1, R = ins;
    int i_lo = ins, i_hi = ins - 1;   // empty sentinel (i_hi < i_lo)
    for (int k = 0; k < cap; ++k) {
        bool haveL = (L >= lo), haveR = (R <= hi);
        if (!haveL && !haveR) break;
        bool takeL = haveL && (!haveR || (center - bp[L]) <= (bp[R] - center));
        if (takeL) { i_lo = L; --L; } else { i_hi = R; ++R; }
    }
    if (i_hi < i_lo) return std::make_pair(0, 0);
    return std::make_pair(i_lo, i_hi + 1);
}

// Keep at most `cap` columns of a Cell by taking an EVEN STRIDE through it.
// Thinning (rather than truncating to the first `cap`) preserves coverage of
// the whole bp window, so the background GRM still represents the region.
static void thin_cell(Cell& c, int cap) {
    const int m = (int) c.X.cols();
    if (cap <= 0 || m <= cap) return;
    const int stride = (m + cap - 1) / cap;
    std::vector<int> keep;
    for (int j = 0; j < m && (int) keep.size() < cap; j += stride) keep.push_back(j);
    GenoMat Z(c.X.rows(), (int) keep.size());
    std::vector<float> mf(keep.size());
    for (size_t k = 0; k < keep.size(); ++k) { Z.col((int) k) = c.X.col(keep[k]); mf[k] = c.maf[keep[k]]; }
    c.X = std::move(Z); c.maf = mf;
}

// Assemble the region for one chunk: [target | flank(both sides) | common].
static bool make_chunk(const ChunkContext& ctx, size_t ci, int a, int b,
                       int chr_lo, int chr_hi, int cc_lo, int cc_hi,
                       int flank_snps, long common_bp, bool common_window_given,
                       int max_common, ChunkData& cd) {
    cd = ChunkData();
    cd.chr = ctx.chr_order[ci];
    cd.start = ctx.wes_bim.bp[a]; cd.end = ctx.wes_bim.bp[b - 1];
    cd.i0_target = a;

    const int fL0 = std::max(chr_lo, a - flank_snps), fL1 = a;
    const int fR0 = b, fR1 = std::min(chr_hi + 1, b + flank_snps);

    const std::vector<float>* wp = ctx.snp_weights.empty() ? 0 : &ctx.snp_weights;
    Cell tgt;
    if (!read_cell_idx(ctx.wes, ctx.geno_keep, a, b, tgt)) return false;
    if (ctx.covZ.cols() > 0) project_covariates(tgt.X, ctx.covZ, ctx.covM);
    apply_alpha(tgt.X, tgt.maf, ctx.alpha, wp, a);
    cd.m_t = (int) tgt.X.cols();
    if (cd.m_t <= 0) return false;

    Cell fl; fl.X = GenoMat(ctx.n_inds, 0);
    Cell fr; fr.X = GenoMat(ctx.n_inds, 0);
    if (fL1 > fL0) { if (!read_cell_idx(ctx.wes, ctx.geno_keep, fL0, fL1, fl)) return false;
                     if (ctx.covZ.cols() > 0) project_covariates(fl.X, ctx.covZ, ctx.covM);
                     apply_alpha(fl.X, fl.maf, ctx.alpha, wp, fL0); }
    if (fR1 > fR0) { if (!read_cell_idx(ctx.wes, ctx.geno_keep, fR0, fR1, fr)) return false;
                     if (ctx.covZ.cols() > 0) project_covariates(fr.X, ctx.covZ, ctx.covM);
                     apply_alpha(fr.X, fr.maf, ctx.alpha, wp, fR0); }
    cd.m_f = (int) fl.X.cols() + (int) fr.X.cols();
    // NOTE: the "no background" bail-out now happens AFTER the common set is read,
    // because with flank_chunks = 0 the common SNPs take over as the background.

    Cell com; com.X = GenoMat(ctx.n_inds, 0);
    if (ctx.use_common) {
        long ctr = (cd.start + cd.end) / 2;
        if (common_window_given) {
            // FIXED bp window centred on the TARGET chunk -- NOT the span of the
            // chunk+flanks. SNP spacing in exome data is wildly uneven (consecutive
            // SNPs in this dataset can be ~800 kb apart), so the span of a fixed
            // NUMBER of SNPs is unbounded: a region straddling a couple of
            // intergenic gaps would drag in thousands of common SNPs and push K_tot
            // straight back to where the bp-window designs failed. A fixed window
            // plus a hard cap keeps K_tot predictable, which is the entire point of
            // the chunk design.
            long half = common_bp / 2;
            long lo_bp = (ctr > half) ? (ctr - half) : 0;
            if (!read_cell(ctx.com, ctx.common_bim, ctx.common_keep, cc_lo, cc_hi,
                           lo_bp, common_bp, com)) return false;
            thin_cell(com, max_common);
        } else {
            // No window specified: take the max_common nearest common SNPs to the
            // chunk midpoint directly. K_tot is bounded by max_common exactly, so
            // no thinning step is needed.
            std::pair<int,int> rng = nearest_snp_range(ctx.common_bim.bp, cc_lo, cc_hi, ctr, max_common);
            if (rng.second > rng.first)
                if (!read_cell_idx(ctx.com, ctx.common_keep, rng.first, rng.second, com)) return false;
        }
        if (com.X.cols() > 0) {
            if (ctx.covZ.cols() > 0) project_covariates(com.X, ctx.covZ, ctx.covM);
            apply_alpha(com.X, com.maf, ctx.alpha_common, 0, 0);
        }
    }
    cd.m_c = (int) com.X.cols();

    // ---- flank_chunks = 0: the common SNPs BECOME the background component ----
    // The model is then { target, common, sigma_e I } rather than
    // { target, flank, common, sigma_e I }. For rare variants this costs nothing:
    // cross-chunk LD between rare variants is negligible (mean |r| ~ 0.013 in
    // simulation, flat with distance), so the flank absorbs no signal the common
    // background does not. It halves K, and the Gram is O(n K^2), so it is close to
    // a 4x saving on the dominant cost. Do NOT use this for common-variant targets,
    // where local LD is real and the flank is doing genuine work.
    if (cd.m_f <= 0 && cd.m_c > 0) {
        cd.flank_is_common = true;
        cd.m_f = cd.m_c; cd.m_c = 0;                   // common occupies the flank block
    }
    // If both are zero the model is just { target, sigma_e I }: unconditioned, but
    // still estimable, so the chunk is kept rather than skipped.
    cd.K = cd.m_t + cd.m_f + cd.m_c;

    cd.V = GenoMat(ctx.n_inds, cd.K);
    cd.V.leftCols(cd.m_t) = tgt.X;
    int off = cd.m_t;
    if (cd.flank_is_common) { cd.V.middleCols(off, cd.m_f) = com.X; off += cd.m_f; }
    else {
        if (fl.X.cols() > 0) { cd.V.middleCols(off, fl.X.cols()) = fl.X; off += (int) fl.X.cols(); }
        if (fr.X.cols() > 0) { cd.V.middleCols(off, fr.X.cols()) = fr.X; off += (int) fr.X.cols(); }
        if (cd.m_c > 0) cd.V.rightCols(cd.m_c) = com.X;
    }

    cd.cj.resize(cd.K);
    for (int j = 0; j < cd.K; ++j) cd.cj[j] = cd.V.col(j).cast<double>().squaredNorm();
    return true;
}

// A unit of work: one chunk, identified by its SNP-index range on a chromosome.
struct Job { size_t ci; int a, b, chr_lo, chr_hi; };

// ===========================================================================
// PARALLELIZATION
//
// Previously the driver built a batch of ChunkData SERIALLY on the main thread
// and then ran only test_chunk() in parallel. That had three costs:
//   * the whole read/decode/standardize stage (measured 0.39 s per chunk at
//     n = 1e5, i.e. ~1.3 h over 11.7k chunks) was a hard Amdahl floor that no
//     thread count could touch;
//   * threads idled at a barrier during every assembly phase; and
//   * peak memory was batch_size x sizeof(V) = 16 x 467 MB = 7.5 GB, because the
//     whole batch was materialized before any of it was consumed.
//
// Now each worker assembles AND tests its own chunk, so assembly is parallel,
// there is no mid-batch barrier, and only n_threads copies of V are ever live at
// once (batch_size now controls output-flush granularity, not memory). This is
// only possible because BedReader replaced readBedBlock's Rcpp allocation --
// nothing in the worker touches the R API.
// ===========================================================================
struct ChunkWorker : public RcppParallel::Worker {
    const ChunkContext& ctx; const std::vector<Job>& jobs; size_t job0;
    const Eigen::MatrixXd& Y; const GenoMat& Yf;
    const std::vector<double>& Vp; const std::vector<double>& yty;
    const ChunkParams& pr; int flank_snps;
    std::vector<ChunkResult>& out; std::vector<char>& ok;
    ChunkWorker(const ChunkContext& ctx, const std::vector<Job>& jobs, size_t job0,
                const Eigen::MatrixXd& Y, const GenoMat& Yf,
                const std::vector<double>& Vp, const std::vector<double>& yty,
                const ChunkParams& pr, int flank_snps,
                std::vector<ChunkResult>& out, std::vector<char>& ok)
        : ctx(ctx), jobs(jobs), job0(job0), Y(Y), Yf(Yf), Vp(Vp), yty(yty),
          pr(pr), flank_snps(flank_snps), out(out), ok(ok) {}
    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t w = begin; w < end; ++w) {
            const Job& j = jobs[job0 + w];
            int cc_lo, cc_hi; common_chr_range(ctx, ctx.chr_order[j.ci], cc_lo, cc_hi);
            ChunkData cd;
            if (!make_chunk(ctx, j.ci, j.a, j.b, j.chr_lo, j.chr_hi, cc_lo, cc_hi,
                            flank_snps, pr.common_bp, pr.common_window_given,
                            pr.max_common_snps, cd)) { ok[w] = 0; continue; }
            test_chunk(cd, Y, Yf, Vp, yty, pr.spa, pr.spa_thresh, out[w]);
            ok[w] = 1;
        }
    }
};

// Peak memory is now n_threads x (V + the K x K working set), NOT
// batch_size x V as before, because each worker builds and discards its own
// chunk. Report it so the thread count can be chosen against available RAM.
static void report_parallel_plan(const ChunkContext& ctx, const ChunkParams& pr,
                                 int Kmax, int n_comp) {
    int th = pr.n_threads;
    if (th <= 0) {
        const char* e = std::getenv("RCPP_PARALLEL_NUM_THREADS");
        th = (e != 0) ? std::atoi(e) : 0;
        if (th <= 0) th = 1;                      // unknown; report per-thread cost
    }
    const double V_mb  = (double) ctx.n_inds * Kmax * 4.0 / 1048576.0;
    const double KK_mb = (double) Kmax * Kmax * 8.0 / 1048576.0;
    const double per   = V_mb + KK_mb * (n_comp + 4);   // V + W cache + S,A,N,Asym
    Rcout << "  Parallel plan: assembly AND testing now run inside the workers.\n"
          << "  Approx peak memory ~ " << (int) per << " MB per thread"
          << " (V " << (int) V_mb << " MB + " << (int) (KK_mb * (n_comp + 4)) << " MB working set).\n";
    if (pr.n_threads > 0) {
        Rcout << "  With n_threads = " << th << " that is ~" << (per * th / 1024.0)
              << " GB total.\n";
    } else {
        // IMPORTANT: memory is now bounded by the THREAD COUNT, not by batch_size
        // as it was before. Left unset, RcppParallel uses every core, so on a
        // many-core node this can be a large multiple of the per-thread figure.
        Rcout << "  WARNING: n_threads is unset, so all available cores will be used and\n"
              << "  peak memory scales with the core count (~" << (per * 16.0 / 1024.0)
              << " GB at 16 threads, ~" << (per * 64.0 / 1024.0) << " GB at 64).\n"
              << "  Set n_threads explicitly to bound it.\n";
    }
}

static Rcpp::List chunk_driver(ChunkContext& ctx, const ChunkParams& pr) {
    const int P = ctx.n_pheno; const Eigen::MatrixXd& Y = ctx.Y;
    Eigen::setNbThreads(1);
    if (pr.n_threads > 0) {
        std::string nt = std::to_string(pr.n_threads);
#ifdef _WIN32
        _putenv_s("RCPP_PARALLEL_NUM_THREADS", nt.c_str());
#else
        setenv("RCPP_PARALLEL_NUM_THREADS", nt.c_str(), 1);
#endif
    }
    std::vector<double> Vp(P), yty(P);
    for (int t = 0; t < P; ++t) { double m = Y.col(t).mean();
        Vp[t]  = (Y.col(t).array() - m).square().sum() / (Y.rows() - 1);
        yty[t] = Y.col(t).squaredNorm(); }
    const GenoMat Yf = Y.cast<float>();      // cast once, not per trait per chunk

    std::ofstream fout; bool tofile = !pr.out_file.empty();
    auto wr = [](std::ofstream& f, double v) { if (std::isnan(v)) f << "NA"; else f << v; };
    if (tofile) {
        fout.open(pr.out_file.c_str());
        if (!fout.is_open()) stop("Could not open out_file: " + pr.out_file);
        fout << "chr\tstart\tend\tm_chunk\tm_flank\tm_common\tphenotype\tvg\tse_vg\th2";
        if (pr.spa) fout << "\tp_spa\tspa_used";
        fout << "\n";
    }

    // enumerate chunks: consecutive blocks of chunk_size SNPs within a chromosome
    std::vector<Job> jobs;
    for (size_t ci = 0; ci < ctx.chr_order.size(); ++ci) {
        int lo = ctx.chr_lo[ci], hi = ctx.chr_hi[ci];
        for (int a = lo; a <= hi; a += pr.chunk_size) {
            int b = std::min(hi + 1, a + pr.chunk_size);
            if (b - a < pr.min_chunk_snps) continue;
            Job j; j.ci = ci; j.a = a; j.b = b; j.chr_lo = lo; j.chr_hi = hi;
            jobs.push_back(j);
        }
    }
    const int flank_snps = pr.chunk_size * pr.flank_chunks;
    Rcout << "Chunks to test: " << jobs.size()
          << "   (chunk_size " << pr.chunk_size << ", flank " << flank_snps << " SNPs/side)\n";
    const int Kmax = pr.chunk_size * (1 + 2 * pr.flank_chunks) + pr.max_common_snps;
    Rcout << "  K_tot per chunk <= " << Kmax
          << "  (" << (pr.chunk_size * (1 + 2 * pr.flank_chunks)) << " WES + <= "
          << pr.max_common_snps << " common)\n";
    if (pr.flank_chunks == 0)
        Rcout << "  flank_chunks = 0: the common SNPs are the background component.\n"
              << "  Model is { target, common, sigma_e I }. Since the Gram is O(n K^2),\n"
              << "  this is ~" << (double)(pr.chunk_size*3 + pr.max_common_snps)
                                  *(pr.chunk_size*3 + pr.max_common_snps)
                                  /((double)Kmax*Kmax)
              << "x cheaper than flank_chunks = 1. Appropriate for RARE targets only.\n";
    report_parallel_plan(ctx, pr, Kmax, 3);
    if (ctx.use_common) {
        if (pr.common_window_given)
            Rcout << "  Common GRM: fixed " << pr.common_bp
                  << " bp window centred on each chunk, thinned to <= "
                  << pr.max_common_snps << " SNPs -- so K_tot is BOUNDED regardless\n"
                  << "  of how far the chunk's SNPs happen to span.\n";
        else
            Rcout << "  Common GRM: nearest " << pr.max_common_snps
                  << " common SNPs to each chunk's midpoint (no bp window, no thinning).\n";
    }

    long n_done = 0, n_skip = 0; size_t idx = 0;
    typedef std::chrono::steady_clock clk;
    clk::time_point t0 = clk::now();
    while (idx < jobs.size()) {
        // Assembly now happens INSIDE the workers, so a "batch" is just the span
        // of jobs handed to one parallelFor -- it controls output flushing and
        // progress granularity, not peak memory.
        const size_t nb = std::min((size_t) std::max(1, pr.batch_size), jobs.size() - idx);
        std::vector<ChunkResult> res(nb);
        std::vector<char> ok(nb, 0);
        ChunkWorker wk(ctx, jobs, idx, Y, Yf, Vp, yty, pr, flank_snps, res, ok);
        RcppParallel::parallelFor(0, nb, wk);
        idx += nb;

        for (size_t b = 0; b < res.size(); ++b) {
            const ChunkResult& cr = res[b];
            if (!ok[b]) { ++n_skip; continue; }
            ++n_done;
            if (!tofile || cr.vg.empty()) continue;
            for (int t = 0; t < P; ++t) {
                fout << cr.chr << '\t' << cr.start << '\t' << cr.end << '\t'
                     << cr.m_t << '\t' << cr.m_f << '\t' << cr.m_c << '\t'
                     << as<std::string>(ctx.trait_names[t]) << '\t';
                wr(fout, cr.vg[t]);    fout << '\t';
                wr(fout, cr.se_vg[t]); fout << '\t';
                wr(fout, cr.h2[t]);
                if (pr.spa) { fout << '\t'; wr(fout, cr.p_spa[t]); fout << '\t' << cr.spa_used[t]; }
                fout << '\n';
            }
        }
        if (tofile) fout.flush();
        double el = std::chrono::duration<double>(clk::now() - t0).count();
        double seen = (double)(n_done + n_skip);
        double eta = (seen > 0) ? el * ((double) jobs.size() - seen) / seen : 0.0;
        Rcout << "\r[chunk] " << (long) seen << "/" << jobs.size()
              << "  " << (int)(100.0 * seen / jobs.size()) << "%"
              << "  elapsed " << (int) el << "s  eta " << (int) eta << "s     " << std::flush;
        Rcpp::checkUserInterrupt();
    }
    Rcout << "\n";
    if (tofile) fout.close();
    Rcout << "Chunks tested: " << n_done << "\n";
    if (n_skip > 0) Rcout << "Chunks skipped (too few SNPs, no background, or read error): " << n_skip << "\n";
    return List::create(_["n_chunks"] = (double) n_done, _["trait_names"] = ctx.trait_names);
}

// ===========================================================================
// ANNOTATION PATH (only used when an annotation matrix is supplied).
//
// Instead of one target = the whole chunk, the chunk's middle SNPs are split
// by the annotation into functional CATEGORIES -- EVERY annotation column is a
// tested category. The background is NOT taken from the annotation: it is built
// purely from the adjacent chunks (flank_chunks) and, if supplied, the common
// SNPs. The model per chunk is:
//   { cat_0, cat_1, ..., cat_{A-1}  (each SPA-tested),
//     flank (neighbouring chunks + any middle SNP in no category),
//     [common], sigma_e I }
// With flank_chunks = 0 and a complete annotation the flank block is empty, and
// the common SNPs become the background component instead.
// and one SPA p-value is produced per category. The category matrices are laid
// out as contiguous leading blocks of V so each is a single component, exactly
// as the single target block was. Everything else (trace normalization, the
// Wald pre-screen, the saddlepoint) is identical to test_chunk, just looped
// over the categories.
// ===========================================================================
struct ChunkDataA {
    std::string chr; long start, end;
    std::vector<int> cat_m;            // columns per tested category (V leading blocks)
    std::vector<std::string> cat_name;
    int m_flank, m_common, K;
    GenoMat V;                         // [cat_0 | ... | cat_{A-1} | flank | common]
    std::vector<double> cj;
    bool flank_is_common;              // flank_chunks = 0: common SNPs are the background
    ChunkDataA() : m_flank(0), m_common(0), K(0), flank_is_common(false) {}
};
struct ChunkResultA {
    std::string chr; long start, end;
    int m_flank, m_common;
    std::vector<std::string> cat_name;
    std::vector<int> cat_m;
    std::vector< std::vector<double> > vg, se_vg, h2, p_spa;   // [cat][trait]
    std::vector< std::vector<int> > spa_used;                   // [cat][trait]
};

// Assemble one chunk, partitioning the target's columns by annotation category.
static bool make_chunk_annot(const ChunkContext& ctx, size_t ci, int a, int b,
                             int chr_lo, int chr_hi, int cc_lo, int cc_hi,
                             int flank_snps, long common_bp, bool common_window_given,
                             int max_common, ChunkDataA& cd) {
    cd = ChunkDataA();
    cd.chr = ctx.chr_order[ci];
    cd.start = ctx.wes_bim.bp[a]; cd.end = ctx.wes_bim.bp[b - 1];
    const std::vector<float>* wp = ctx.snp_weights.empty() ? 0 : &ctx.snp_weights;

    Cell tgt;
    if (!read_cell_idx(ctx.wes, ctx.geno_keep, a, b, tgt)) return false;
    if (ctx.covZ.cols() > 0) project_covariates(tgt.X, ctx.covZ, ctx.covM);
    apply_alpha(tgt.X, tgt.maf, ctx.alpha, wp, a);
    int m_t = (int) tgt.X.cols();
    if (m_t <= 0) return false;

    // partition the chunk's columns by category id (-1 = background -> flank)
    int A = ctx.n_annot_cat;
    std::vector< std::vector<int> > cat_cols(A);
    std::vector<int> bg_cols;
    for (int j = 0; j < m_t; ++j) {
        int cid = (a + j >= 0 && a + j < (int) ctx.snp_cat.size()) ? ctx.snp_cat[a + j] : -1;
        if (cid >= 0 && cid < A) cat_cols[cid].push_back(j);
        else bg_cols.push_back(j);
    }
    std::vector<int> act;
    for (int c = 0; c < A; ++c) if (!cat_cols[c].empty()) act.push_back(c);
    if (act.empty()) return false;                 // no category SNPs in this chunk to test

    // neighbouring flank chunks (same as the no-annotation path)
    const int fL0 = std::max(chr_lo, a - flank_snps), fL1 = a;
    const int fR0 = b, fR1 = std::min(chr_hi + 1, b + flank_snps);
    Cell fl; fl.X = GenoMat(ctx.n_inds, 0);
    Cell fr; fr.X = GenoMat(ctx.n_inds, 0);
    if (fL1 > fL0) { if (!read_cell_idx(ctx.wes, ctx.geno_keep, fL0, fL1, fl)) return false;
                     if (ctx.covZ.cols() > 0) project_covariates(fl.X, ctx.covZ, ctx.covM);
                     apply_alpha(fl.X, fl.maf, ctx.alpha, wp, fL0); }
    if (fR1 > fR0) { if (!read_cell_idx(ctx.wes, ctx.geno_keep, fR0, fR1, fr)) return false;
                     if (ctx.covZ.cols() > 0) project_covariates(fr.X, ctx.covZ, ctx.covM);
                     apply_alpha(fr.X, fr.maf, ctx.alpha, wp, fR0); }
    // Background = adjacent chunks + any middle SNP in no category. The bail-out
    // moved BELOW the common read: with flank_chunks = 0 and a complete annotation
    // this count is legitimately zero, and the common SNPs take over as the
    // background. (Previously that combination made every chunk return false, so
    // the run wrote a header and then silently computed nothing.)
    int m_flank = (int) bg_cols.size() + (int) fl.X.cols() + (int) fr.X.cols();

    // common GRM (fixed bp window centred on the chunk, thinned; or the nearest
    // max_common SNPs directly if no window is given) -- as no-annot
    Cell com; com.X = GenoMat(ctx.n_inds, 0);
    if (ctx.use_common) {
        long ctr = (cd.start + cd.end) / 2;
        if (common_window_given) {
            long half = common_bp / 2;
            long lo_bp = (ctr > half) ? (ctr - half) : 0;
            if (!read_cell(ctx.com, ctx.common_bim, ctx.common_keep, cc_lo, cc_hi,
                           lo_bp, common_bp, com)) return false;
            thin_cell(com, max_common);
        } else {
            std::pair<int,int> rng = nearest_snp_range(ctx.common_bim.bp, cc_lo, cc_hi, ctr, max_common);
            if (rng.second > rng.first)
                if (!read_cell_idx(ctx.com, ctx.common_keep, rng.first, rng.second, com)) return false;
        }
        if (com.X.cols() > 0) {
            if (ctx.covZ.cols() > 0) project_covariates(com.X, ctx.covZ, ctx.covM);
            apply_alpha(com.X, com.maf, ctx.alpha_common, 0, 0);
        }
    }
    int m_c = (int) com.X.cols();

    // flank_chunks = 0 (or no uncategorized middle SNPs): the common set becomes
    // the background component, exactly as in the no-annotation path.
    if (m_flank <= 0 && m_c > 0) {
        cd.flank_is_common = true;
        m_flank = m_c; m_c = 0;
    }
    // If both are zero the model is { cat_0, ..., cat_{A-1}, sigma_e I }: each
    // category is conditioned only on the other categories. Unconditioned on any
    // regional background, but still estimable, so keep the chunk.

    int Kt = 0; for (size_t k = 0; k < act.size(); ++k) Kt += (int) cat_cols[act[k]].size();
    cd.K = Kt + m_flank + m_c;
    cd.V = GenoMat(ctx.n_inds, cd.K);
    int off = 0;
    for (size_t k = 0; k < act.size(); ++k) {
        int c = act[k]; const std::vector<int>& cc = cat_cols[c];
        for (size_t q = 0; q < cc.size(); ++q) cd.V.col(off + (int) q) = tgt.X.col(cc[q]);
        off += (int) cc.size();
        cd.cat_m.push_back((int) cc.size());
        cd.cat_name.push_back(ctx.annot_names[c]);
    }
    if (cd.flank_is_common) {
        cd.V.middleCols(off, m_flank) = com.X; off += m_flank;
    } else {
        for (size_t q = 0; q < bg_cols.size(); ++q) cd.V.col(off++) = tgt.X.col(bg_cols[q]);
        if (fl.X.cols() > 0) { cd.V.middleCols(off, fl.X.cols()) = fl.X; off += (int) fl.X.cols(); }
        if (fr.X.cols() > 0) { cd.V.middleCols(off, fr.X.cols()) = fr.X; off += (int) fr.X.cols(); }
        if (m_c > 0) cd.V.rightCols(m_c) = com.X;
    }
    cd.m_flank = m_flank; cd.m_common = m_c;
    cd.cj.resize(cd.K);
    for (int j = 0; j < cd.K; ++j) cd.cj[j] = cd.V.col(j).cast<double>().squaredNorm();
    return true;
}

// Multi-category SPA test. Components: 0..A-1 categories (tested), A flank,
// then flank (if any), then common (if any), env residual. Mirrors test_chunk but loops the SPA over
// the A categories, reusing the single per-chunk Cholesky.
static void test_chunk_annot(const ChunkDataA& cd, const Eigen::MatrixXd& Y, const GenoMat& Yf,
                             const std::vector<double>& Vp, const std::vector<double>& yty,
                             bool spa, double spa_thresh,
                             ChunkResultA& cr) {
    const int n = (int) Y.rows(), P = (int) Y.cols(), K = cd.K;
    const int A = (int) cd.cat_m.size();
    const bool has_f = (cd.m_flank  > 0);
    const bool has_c = (cd.m_common > 0);
    // Components: the A tested categories, then flank and common ONLY if present.
    // Built this way so flank_chunks = 0 with a complete annotation and no common
    // fileset still runs, as { cat_0, ..., cat_{A-1}, sigma_e I }.
    const int C = A + (has_f ? 1 : 0) + (has_c ? 1 : 0), env = C;

    cr.chr = cd.chr; cr.start = cd.start; cr.end = cd.end;
    // Report honestly when the common set is standing in as the background.
    cr.m_flank  = cd.flank_is_common ? 0           : cd.m_flank;
    cr.m_common = cd.flank_is_common ? cd.m_flank  : cd.m_common;
    cr.cat_name = cd.cat_name; cr.cat_m = cd.cat_m;
    cr.vg.assign(A, std::vector<double>(P, NA_REAL));
    cr.se_vg.assign(A, std::vector<double>(P, NA_REAL));
    cr.h2.assign(A, std::vector<double>(P, NA_REAL));
    cr.p_spa.assign(A, std::vector<double>(P, NA_REAL));
    cr.spa_used.assign(A, std::vector<int>(P, 0));
    if (A <= 0) return;

    // component column offsets in V: cats, then flank (if any), then common (if any)
    std::vector<int> off(C + 1, 0);
    for (int c = 0; c < A; ++c) off[c + 1] = off[c] + cd.cat_m[c];
    int nx = A;
    if (has_f) { off[nx + 1] = off[nx] + cd.m_flank;  ++nx; }
    if (has_c) { off[nx + 1] = off[nx] + cd.m_common; ++nx; }

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
    for (int a = 0; a < C; ++a) { T(a, env) = (double) n; T(env, a) = (double) n; }
    T(env, env) = (double) n;
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

    Eigen::MatrixXd U = (cd.V.transpose() * Yf).cast<double>();      // K x P, one GEMM

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
            cr.se_vg[c][t] = std::sqrt(var0);
            double p_wald = std::erfc(std::abs(sigma[c] / cr.se_vg[c][t]) / std::sqrt(2.0));
            if (p_wald >= spa_thresh) { cr.p_spa[c][t] = p_wald; cr.spa_used[c][t] = 0; continue; }

            Eigen::LLT<Eigen::MatrixXd> lltS(Smat);
            if (lltS.info() != Eigen::Success) { cr.p_spa[c][t] = p_wald; continue; }
            Eigen::MatrixXd Cf = lltS.matrixL();
            Eigen::MatrixXd Asym = Cf.transpose() * Amat * Cf;
            Asym = (0.5 * (Asym + Asym.transpose())).eval();
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ses(Asym, Eigen::EigenvaluesOnly);
            if (ses.info() != Eigen::Success) { cr.p_spa[c][t] = p_wald; continue; }
            const Eigen::VectorXd& ev = ses.eigenvalues();
            std::vector<double> eig(ev.data(), ev.data() + K);
            QuadSpaResult qr = quad_spa_solve(sigma[c], eig, eig_rep, n_rep);
            if (qr.converged) { cr.p_spa[c][t] = qr.p; cr.spa_used[c][t] = 1; }
            else              { cr.p_spa[c][t] = p_wald; cr.spa_used[c][t] = 0; }
        }
    }
}

struct ChunkWorkerA : public RcppParallel::Worker {
    const ChunkContext& ctx; const std::vector<Job>& jobs; size_t job0;
    const Eigen::MatrixXd& Y; const GenoMat& Yf;
    const std::vector<double>& Vp; const std::vector<double>& yty;
    const ChunkParams& pr; int flank_snps;
    std::vector<ChunkResultA>& out; std::vector<char>& ok;
    ChunkWorkerA(const ChunkContext& ctx, const std::vector<Job>& jobs, size_t job0,
                 const Eigen::MatrixXd& Y, const GenoMat& Yf,
                 const std::vector<double>& Vp, const std::vector<double>& yty,
                 const ChunkParams& pr, int flank_snps,
                 std::vector<ChunkResultA>& out, std::vector<char>& ok)
        : ctx(ctx), jobs(jobs), job0(job0), Y(Y), Yf(Yf), Vp(Vp), yty(yty),
          pr(pr), flank_snps(flank_snps), out(out), ok(ok) {}
    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t w = begin; w < end; ++w) {
            const Job& j = jobs[job0 + w];
            int cc_lo, cc_hi; common_chr_range(ctx, ctx.chr_order[j.ci], cc_lo, cc_hi);
            ChunkDataA cd;
            if (!make_chunk_annot(ctx, j.ci, j.a, j.b, j.chr_lo, j.chr_hi, cc_lo, cc_hi,
                                  flank_snps, pr.common_bp, pr.common_window_given,
                                  pr.max_common_snps, cd)) { ok[w] = 0; continue; }
            test_chunk_annot(cd, Y, Yf, Vp, yty, pr.spa, pr.spa_thresh, out[w]);
            ok[w] = 1;
        }
    }
};

static Rcpp::List chunk_driver_annot(ChunkContext& ctx, const ChunkParams& pr) {
    const int P = ctx.n_pheno; const Eigen::MatrixXd& Y = ctx.Y;
    Eigen::setNbThreads(1);
    if (pr.n_threads > 0) {
        std::string nt = std::to_string(pr.n_threads);
#ifdef _WIN32
        _putenv_s("RCPP_PARALLEL_NUM_THREADS", nt.c_str());
#else
        setenv("RCPP_PARALLEL_NUM_THREADS", nt.c_str(), 1);
#endif
    }
    std::vector<double> Vp(P), yty(P);
    for (int t = 0; t < P; ++t) { double m = Y.col(t).mean();
        Vp[t]  = (Y.col(t).array() - m).square().sum() / (Y.rows() - 1);
        yty[t] = Y.col(t).squaredNorm(); }
    const GenoMat Yf = Y.cast<float>();

    std::ofstream fout; bool tofile = !pr.out_file.empty();
    auto wr = [](std::ofstream& f, double v) { if (std::isnan(v)) f << "NA"; else f << v; };
    if (tofile) {
        fout.open(pr.out_file.c_str());
        if (!fout.is_open()) stop("Could not open out_file: " + pr.out_file);
        fout << "chr\tstart\tend\tcategory\tm_cat\tm_flank\tm_common\tphenotype\tvg\tse_vg\th2";
        if (pr.spa) fout << "\tp_spa\tspa_used";
        fout << "\n";
    }

    std::vector<Job> jobs;
    for (size_t ci = 0; ci < ctx.chr_order.size(); ++ci) {
        int lo = ctx.chr_lo[ci], hi = ctx.chr_hi[ci];
        for (int a = lo; a <= hi; a += pr.chunk_size) {
            int b = std::min(hi + 1, a + pr.chunk_size);
            if (b - a < pr.min_chunk_snps) continue;
            Job j; j.ci = ci; j.a = a; j.b = b; j.chr_lo = lo; j.chr_hi = hi;
            jobs.push_back(j);
        }
    }
    const int flank_snps = pr.chunk_size * pr.flank_chunks;
    Rcout << "Chunks to test: " << jobs.size() << "   (chunk_size " << pr.chunk_size
          << ", flank " << flank_snps << " SNPs/side, " << ctx.n_annot_cat << " categories)\n";
    report_parallel_plan(ctx, pr, pr.chunk_size * (1 + 2 * pr.flank_chunks) + pr.max_common_snps,
                         ctx.n_annot_cat + 2);

    long n_done = 0, n_skip = 0; size_t idx = 0;
    typedef std::chrono::steady_clock clk;
    clk::time_point t0 = clk::now();
    while (idx < jobs.size()) {
        const size_t nb = std::min((size_t) std::max(1, pr.batch_size), jobs.size() - idx);
        std::vector<ChunkResultA> res(nb);
        std::vector<char> ok(nb, 0);
        ChunkWorkerA wk(ctx, jobs, idx, Y, Yf, Vp, yty, pr, flank_snps, res, ok);
        RcppParallel::parallelFor(0, nb, wk);
        idx += nb;

        for (size_t bb = 0; bb < res.size(); ++bb) {
            const ChunkResultA& cr = res[bb];
            if (!ok[bb]) { ++n_skip; continue; }
            ++n_done;
            if (!tofile || cr.cat_name.empty()) continue;
            int A = (int) cr.cat_name.size();
            for (int c = 0; c < A; ++c)
                for (int t = 0; t < P; ++t) {
                    fout << cr.chr << '\t' << cr.start << '\t' << cr.end << '\t'
                         << cr.cat_name[c] << '\t' << cr.cat_m[c] << '\t'
                         << cr.m_flank << '\t' << cr.m_common << '\t'
                         << as<std::string>(ctx.trait_names[t]) << '\t';
                    wr(fout, cr.vg[c][t]);    fout << '\t';
                    wr(fout, cr.se_vg[c][t]); fout << '\t';
                    wr(fout, cr.h2[c][t]);
                    if (pr.spa) { fout << '\t'; wr(fout, cr.p_spa[c][t]); fout << '\t' << cr.spa_used[c][t]; }
                    fout << '\n';
                }
        }
        if (tofile) fout.flush();
        double el = std::chrono::duration<double>(clk::now() - t0).count();
        double seen = (double)(n_done + n_skip);
        double eta = (seen > 0) ? el * ((double) jobs.size() - seen) / seen : 0.0;
        Rcout << "\r[chunk-annot] " << (long) seen << "/" << jobs.size()
              << "  " << (int)(100.0 * seen / jobs.size()) << "%"
              << "  elapsed " << (int) el << "s  eta " << (int) eta << "s     " << std::flush;
        Rcpp::checkUserInterrupt();
    }
    Rcout << "\n";
    if (tofile) fout.close();
    Rcout << "Chunks tested: " << n_done << "\n";
    if (n_skip > 0) Rcout << "Chunks skipped (no category SNPs, no background, or read error): " << n_skip << "\n";
    return List::create(_["n_chunks"] = (double) n_done, _["trait_names"] = ctx.trait_names,
                        _["categories"] = wrap(ctx.annot_names));
}

}  // end anonymous namespace

// ===========================================================================
// Export
// ===========================================================================
// [[Rcpp::export]]
Rcpp::List he_chunk_spa(const std::string& filename,
                        const SEXP pheno_mat,
                        int chunk_size = 256,
                        int flank_chunks = 1,
                        int min_chunk_snps = 5,
                        double alpha = -1.0,
                        double alpha_common = -1.0,
                        Rcpp::Nullable<Rcpp::String> common_filename = R_NilValue,
                        double common_window = NA_REAL,
                        int max_common_snps = 256,
                        std::string out_file = "",
                        int batch_size = 64,
                        int n_threads = 0,
                        Rcpp::Nullable<Rcpp::NumericVector> weights = R_NilValue,
                        Rcpp::Nullable<Rcpp::NumericMatrix> covariates = R_NilValue,
                        bool SPA = true,
                        double spa_pval_threshold = 0.1,
                        Rcpp::Nullable<Rcpp::IntegerMatrix> annotation = R_NilValue,
                        Rcpp::Nullable<Rcpp::CharacterVector> annot_names = R_NilValue) {
    if (chunk_size < 1) stop("chunk_size must be >= 1");
    if (flank_chunks < 0) stop("flank_chunks must be >= 0");
    ChunkContext ctx = setup_chunk_context(filename, pheno_mat, alpha, alpha_common,
                                           common_filename, weights, covariates,
                                           annotation, annot_names);
    ChunkParams pr;
    pr.chunk_size = chunk_size; pr.flank_chunks = flank_chunks;
    pr.min_chunk_snps = min_chunk_snps;
    // common_window left unspecified (NA) -> take the max_common_snps nearest
    // common SNPs to each chunk's midpoint directly, no bp window, no thinning.
    pr.common_window_given = !ISNAN(common_window);
    pr.common_bp = pr.common_window_given ? (long) common_window : 0;
    pr.max_common_snps = max_common_snps;
    pr.spa = SPA; pr.spa_thresh = spa_pval_threshold;
    pr.out_file = out_file; pr.batch_size = batch_size; pr.n_threads = n_threads;
    // Annotation supplied -> per-category SPA path. Every annotation column is a
    // tested category; the background comes from the adjacent chunks and/or the
    // common SNPs, so with flank_chunks = 0 a common fileset is required.
    if (ctx.use_annot) {
        if (flank_chunks < 1 && !ctx.use_common)
            Rcpp::warning("flank_chunks = 0 with no common_filename: the only background "
                          "left is middle SNPs belonging to no annotation category. If the "
                          "annotation covers every SNP, each category is conditioned only on "
                          "the other categories -- estimates will absorb regional background "
                          "and are expected to be inflated.");
        return chunk_driver_annot(ctx, pr);
    }
    // No annotation: the model still needs SOME background component, but that can
    // come either from the neighbouring chunks (flank_chunks >= 1) or, when
    // flank_chunks = 0, from the common-SNP set. One of the two must be present.
    if (flank_chunks < 1 && !ctx.use_common)
        Rcpp::warning("flank_chunks = 0 with no common_filename: there is NO background "
                      "component, so the model is just { target, sigma_e I }. Each chunk's "
                      "variance is unconditioned on its neighbours or on regional structure, "
                      "so vg is expected to be inflated. Intended for testing only.");
    return chunk_driver(ctx, pr);
}
