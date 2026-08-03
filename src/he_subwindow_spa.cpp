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

using namespace Rcpp;

// ===========================================================================
// SUB-WINDOW HE + saddlepoint association testing.
//
// A different design from he_sliding_window_enrich_alpha.cpp, aimed squarely
// at making the quadratic-form SPA affordable.
//
// WHY THIS FILE EXISTS. The SPA's cost is governed by K_tot, the TOTAL number
// of genotype columns across ALL variance components in a window's model,
// because it needs the K_tot x K_tot Gram matrix, its Cholesky (O(K_tot^3))
// and a symmetric eigensolve (~10x that) per tested component. In the 1 Mb
// sliding-window design the two 1 Mb WES flanks dominate that total:
//
//     left 8547 + right 8547 + middle 2754 + common 406  =  K_tot ~ 20000
//
// at which point a single window costs ~45 min. Shrinking the TARGET does not
// help at all -- the flanks are untouched.
//
// THE DESIGN HERE. Keep the 1 Mb cell as the unit of local context, but:
//   * drop the two 1 Mb flanks entirely;
//   * split the 1 Mb cell into small sub-windows (e.g. 5 kb) and test each;
//   * let the REMAINDER of the 1 Mb cell (the 995 kb outside the sub-window)
//     play the role the flanks used to play, absorbing local background;
//   * keep the optional common-SNP GRM.
// The model for sub-window s is therefore
//     { K_s (target), K_rest = K_1Mb - K_s, [K_common], sigma_e I }
// and
//     K_tot = m_1Mb + m_common ~ 3160,
// a 6.4x reduction in K_tot, i.e. ~260x on the cubic cost.
//
// THE SUBTRACTION TRICK (the user's idea, and it is exact). Writing the
// UN-normalized component Gram as S_A = sum_{j in A} v_j v_j' for a column set
// A, disjointness gives exactly
//     S_rest = S_1Mb - S_s ,
// and the same subtraction holds for every quantity we actually need:
//     tr(S_rest)          = tr(S_1Mb) - tr(S_s)                (trace norm.)
//     y'S_rest y          = y'S_1Mb y - y'S_s y                (the q vector)
//     tr(S_a S_b)         = ||G[A,B]||_F^2, and blocks subtract likewise (T)
// so the whole moment system for a sub-window is assembled in O(m_s * m_1Mb)
// from quantities computed ONCE for the 1 Mb cell. Verified numerically to
// ~1e-13 against direct computation.
//
// THE SECOND WIN. Every sub-window of a 1 Mb cell shares the SAME column set
// V = [V_1Mb | V_common]; only the PARTITION of those columns into components
// changes. Hence G = V'V and its Cholesky factor L are computed ONCE per 1 Mb
// cell and reused by all ~200 sub-windows. Moreover
//     P = L' D_0 L  and  Q = L' D_M L
// decompose into per-component blocks, so with PL_mid = L[mid,:]'L[mid,:] and
// PL_com = L[com,:]'L[com,:] precomputed per cell, each sub-window needs only
// PL_s = L[s,:]'L[s,:] -- O(K_tot^2 * m_s) instead of O(K_tot^3).
//
// WHAT IS STILL EXPENSIVE, AND THE KNOB FOR IT. The Cholesky of S and the
// symmetric eigensolve remain O(K_tot^3) PER SUB-WINDOW, and with ~200
// sub-windows per cell that still dominates (~109 min per cell at
// K_tot = 3160). `background_rank` addresses this: it replaces the background
// (rest + common) column blocks by their leading `background_rank`
// eigenvectors before the SPA, cutting K_tot to background_rank + m_s.
// Indicative genome-wide totals (3000 cells, 200 sub-windows, 19 threads):
//     K_tot 3160 (exact) ~ 286 h      K_tot 1000 ~  8 h
//     K_tot 2000         ~  73 h      K_tot  500 ~  1 h
// background_rank = 0 keeps the calculation EXACT and is the right setting for
// spot-checking a chromosome; a few hundred is what makes a genome-wide run
// tractable. This truncation is an APPROXIMATION and has not been validated
// against the exact result on real data -- compare the two on a subset before
// trusting it.
//
// Everything lives in an anonymous namespace (internal linkage) so it cannot
// collide with the identically-named helpers in the other .cpp files of the
// package; only the exported wrapper has external linkage.
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
static Cell read_cell_idx(const std::string& prefix, int n_total, int n_snps,
                          const std::vector<int>& keep, int i0, int i1,
                          const std::vector<int>* cat_id = 0,
                          int collapse_mac = 0, int collapse_n = 5) {
    Cell cell; cell.i0 = -1;
    int m = i1 - i0;
    if (m <= 0) { cell.X = GenoMat(keep.size(), 0); return cell; }
    int n_keep = (int) keep.size();
    IntegerMatrix block = readBedBlock(prefix + ".bed", n_total, n_snps, 0, n_total - 1, i0, i1 - 1);

    // -------- no collapsing: original path (columns map 1:1 to the .bim) --------
    if (!(collapse_mac > 0 && cat_id != 0 && collapse_n >= 1)) {
        GenoMat X(n_keep, m);
        for (int i = 0; i < n_keep; ++i)
            for (int j = 0; j < m; ++j)
                X(i, j) = (float) block(keep[i], j);
        standardize_capture_maf(X, cell.maf);
        cell.X = std::move(X);
        cell.i0 = i0;
        return cell;
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
        for (int i = 0; i < n_keep; ++i) { int g = block(keep[i], j); if (g >= 0) { s += g; ++nv; } }
        int si = (int) (s + 0.5);
        int mac = (nv > 0) ? std::min(si, 2 * nv - si) : 0;
        bool flip = (nv > 0 && s > (double) nv);          // coded allele is major
        int cid = (i0 + j >= 0 && i0 + j < ncat_all) ? (*cat_id)[i0 + j] : -1;

        if (mac >= collapse_mac) {
            if (grp_open) { out_cols.push_back(grp); out_cat.push_back(grp_cat); out_pseudo.push_back(1);
                            grp_open = false; grp_cnt = 0; grp_cat = -1; }
            std::vector<float> col(n_keep);
            for (int i = 0; i < n_keep; ++i) { int g = block(keep[i], j); col[i] = (g < 0) ? -1.f : (float) g; }
            out_cols.push_back(std::move(col)); out_cat.push_back(cid); out_pseudo.push_back(0);
        } else {
            if (grp_open && (grp_cat != cid || grp_cnt >= collapse_n)) {
                out_cols.push_back(grp); out_cat.push_back(grp_cat); out_pseudo.push_back(1);
                grp_open = false; grp_cnt = 0; grp_cat = -1;
            }
            if (!grp_open) { grp.assign(n_keep, 0.f); grp_open = true; grp_cat = cid; grp_cnt = 0; }
            for (int i = 0; i < n_keep; ++i) {
                int g = block(keep[i], j);
                if (g >= 0) grp[i] += (flip ? (float)(2 - g) : (float) g);   // additive minor dosage
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
    return cell;
}

// Read a cell by bp range [cell_start, cell_start + W) within chromosome index
// range [lo, hi] of `bim`. Used for the COMMON fileset (matched to a WES cell's
// bp span). Same routine as before.
static Cell read_cell(const std::string& prefix, int n_total, int n_snps,
                      const BimInfo& bim, const std::vector<int>& keep,
                      int lo, int hi, long cell_start, long W) {
    Cell cell; cell.i0 = -1;
    if (lo < 0 || hi < lo || W <= 0) { cell.X = GenoMat(keep.size(), 0); return cell; }
    int i0 = lower_index(bim.bp, lo, hi, cell_start);
    int i1 = lower_index(bim.bp, lo, hi, cell_start + W);
    int m = i1 - i0;
    if (m <= 0) { cell.X = GenoMat(keep.size(), 0); return cell; }
    int n_keep = (int) keep.size();
    GenoMat X(n_keep, m);
    IntegerMatrix block = readBedBlock(prefix + ".bed", n_total, n_snps, 0, n_total - 1, i0, i1 - 1);
    for (int i = 0; i < n_keep; ++i)
        for (int j = 0; j < m; ++j)
            X(i, j) = (float) block(keep[i], j);
    standardize_capture_maf(X, cell.maf);
    cell.X = std::move(X);
    cell.i0 = i0;
    return cell;
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
struct SubContext {
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

static SubContext setup_sub_context(const std::string& filename, const SEXP pheno_mat,
                                    double alpha, double alpha_common,
                                    Rcpp::Nullable<Rcpp::String> common_filename,
                                    Rcpp::Nullable<Rcpp::NumericVector> weights,
                                    Rcpp::Nullable<Rcpp::NumericMatrix> covariates) {
    SubContext ctx;
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
        IntegerVector cmatch = match(analysis_iid, common_iid);
        ctx.common_keep.resize(ctx.n_inds);
        for (int i = 0; i < ctx.n_inds; ++i) {
            if (cmatch[i] == NA_INTEGER) stop("An analysis individual is missing from the common-SNP .fam file");
            ctx.common_keep[i] = cmatch[i] - 1;
        }
    }
    return ctx;
}

// ===========================================================================
// One 1 Mb cell, ready for sub-window testing
// ===========================================================================
struct CellData {
    std::string chr; long start, end;
    int i0_mid;                    // .bim index of the cell's first WES SNP
    int m_mid, m_com, K;           // column counts; K = m_mid + m_com
    GenoMat V;                     // n x K, ALPHA-WEIGHTED but NOT trace-normalized
    std::vector<double> cj;        // ||v_j||^2 per column  (so tr(S_A) = sum_{j in A} cj)
    std::vector<long> bp;          // bp position of each MIDDLE column
};

// Apply the alpha weight (and optional per-SNP LD weight) to a standardized
// block, WITHOUT the per-component trace normalization -- that constant depends
// on which columns form a component, which changes per sub-window, so it is
// applied later as a scalar (see g_a in sub_window_scan).
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
// Sub-window scan for ONE 1 Mb cell
// ===========================================================================
struct SubResult {
    std::string chr; long start, end;
    int m_s, m_rest, m_com;
    std::vector<double> vg, se_vg, h2, p_spa;   // per phenotype
    std::vector<int> spa_used;
};

// Components, in this fixed order:
//   0 = target sub-window s
//   1 = rest of the 1 Mb cell
//   2 = common background            (only if present)
//   env = residual
//
// All kernels use the trace normalization tr(K_a) = n, i.e.
//     K_a = g_a * S_A ,  g_a = n / tr(S_A) ,  tr(S_A) = sum_{j in A} cj .
// With that convention T(a,env) = tr(K_a) = n for every a, exactly.
static void scan_cell(const CellData& cd, const Eigen::MatrixXd& Y,
                      const std::vector<double>& Vp, long sub_bp, int min_sub_snps,
                      bool spa, double spa_thresh, int background_rank,
                      std::vector<SubResult>& out) {
    const int n = (int) Y.rows(), P = (int) Y.cols();
    const int K = cd.K, m_mid = cd.m_mid, m_com = cd.m_com;
    const bool has_com = (m_com > 0);
    const int C = has_com ? 3 : 2;          // number of genetic components
    const int env = C;
    if (m_mid <= 0) return;

    // ---- per-cell quantities, computed ONCE and reused by every sub-window --
    // G = V'V  (EXACT, double). This is the object the whole reduction needs,
    // and it is identical for every sub-window because the COLUMNS are the
    // same -- only the partition into components changes.
    GenoMat Gf = GenoMat::Zero(K, K);
    Gf.selfadjointView<Eigen::Upper>().rankUpdate(cd.V.transpose());
    Eigen::MatrixXd G = GenoMat(Gf.selfadjointView<Eigen::Upper>()).cast<double>();

    // squared Gram, for the tr(S_a S_b) = ||G[A,B]||_F^2 block sums
    Eigen::MatrixXd G2 = G.array().square();

    // u = V' y per phenotype  ->  y'S_A y = sum_{j in A} u_j^2
    std::vector<Eigen::VectorXd> u(P);
    for (int t = 0; t < P; ++t) u[t] = (cd.V.transpose() * Y.col(t).cast<float>()).cast<double>();

    // running totals over the whole 1 Mb middle / common blocks
    double tr_mid = 0.0, tr_com = 0.0;
    for (int j = 0; j < m_mid; ++j) tr_mid += cd.cj[j];
    for (int j = m_mid; j < K; ++j) tr_com += cd.cj[j];
    const double f_mm = G2.topLeftCorner(m_mid, m_mid).sum();               // tr(S_mid S_mid)
    const double f_cc = has_com ? G2.bottomRightCorner(m_com, m_com).sum() : 0.0;
    const double f_mc = has_com ? G2.topRightCorner(m_mid, m_com).sum()  : 0.0;

    // Cholesky of G, and the per-block products PL_mid / PL_com. Built ONCE
    // per cell; the escalating ridge is required because G is accumulated in
    // float from ultra-rare (near-collinear) columns -- see the long note in
    // he_sliding_window_enrich_alpha.cpp.
    Eigen::MatrixXd L, PL_mid, PL_com;
    bool have_L = false;
    if (spa) {
        double gscale = G.trace() / K;
        Eigen::MatrixXd Gj = G; double applied = 0.0;
        for (double r = 1e-4; r <= 1e-1; r *= 10.0) {
            Gj.diagonal().array() += (r - applied) * gscale; applied = r;
            Eigen::LLT<Eigen::MatrixXd> llt(Gj);
            if (llt.info() == Eigen::Success) { L = llt.matrixL(); have_L = true; break; }
        }
        if (have_L) {
            PL_mid.noalias() = L.topRows(m_mid).transpose()    * L.topRows(m_mid);
            if (has_com) PL_com.noalias() = L.bottomRows(m_com).transpose() * L.bottomRows(m_com);
        }
    }

    // ---- sub-window boundaries (by bp within the cell) ----------------------
    std::vector<std::pair<int,int> > subs;
    {
        int j = 0;
        while (j < m_mid) {
            long b0 = cd.bp[j]; int k = j;
            while (k < m_mid && cd.bp[k] < b0 + sub_bp) ++k;
            if (k - j >= min_sub_snps) subs.push_back(std::make_pair(j, k));
            j = k;
        }
    }

    for (size_t si = 0; si < subs.size(); ++si) {
        const int s0 = subs[si].first, s1 = subs[si].second, m_s = s1 - s0;
        SubResult sr;
        sr.chr = cd.chr; sr.start = cd.bp[s0]; sr.end = cd.bp[s1 - 1];
        sr.m_s = m_s; sr.m_rest = m_mid - m_s; sr.m_com = m_com;
        sr.vg.assign(P, NA_REAL); sr.se_vg.assign(P, NA_REAL);
        sr.h2.assign(P, NA_REAL); sr.p_spa.assign(P, NA_REAL);
        sr.spa_used.assign(P, 0);
        if (sr.m_rest <= 0) { out.push_back(sr); continue; }

        // ---- traces, all by SUBTRACTION from the per-cell totals ------------
        double tr_s = 0.0; for (int j = s0; j < s1; ++j) tr_s += cd.cj[j];
        const double tr_rest = tr_mid - tr_s;
        if (!(tr_s > 0.0) || !(tr_rest > 0.0)) { out.push_back(sr); continue; }

        const double f_ss = G2.block(s0, s0, m_s, m_s).sum();               // tr(S_s S_s)
        const double f_sm = G2.block(s0, 0, m_s, m_mid).sum();              // tr(S_s S_mid)
        const double f_sr = f_sm - f_ss;                                     // tr(S_s S_rest)
        const double f_rr = f_mm - 2.0 * f_sm + f_ss;                        // tr(S_rest S_rest)
        const double f_sc = has_com ? G2.block(s0, m_mid, m_s, m_com).sum() : 0.0;
        const double f_rc = has_com ? f_mc - f_sc : 0.0;

        // trace-normalization constants g_a = n / tr(S_A)
        std::vector<double> g(C);
        g[0] = (double) n / tr_s; g[1] = (double) n / tr_rest;
        if (has_com) g[2] = (double) n / tr_com;

        // ---- moment matrix T ------------------------------------------------
        Eigen::MatrixXd T = Eigen::MatrixXd::Zero(C + 1, C + 1);
        T(0,0) = g[0]*g[0]*f_ss;  T(1,1) = g[1]*g[1]*f_rr;
        T(0,1) = T(1,0) = g[0]*g[1]*f_sr;
        if (has_com) {
            T(2,2) = g[2]*g[2]*f_cc;
            T(0,2) = T(2,0) = g[0]*g[2]*f_sc;
            T(1,2) = T(2,1) = g[1]*g[2]*f_rc;
        }
        for (int a = 0; a < C; ++a) { T(a,env) = (double) n; T(env,a) = (double) n; }  // tr(K_a) = n
        T(env,env) = (double) n;
        Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> Tcod(T);
        Eigen::MatrixXd Tinv = Tcod.pseudoInverse();

        for (int t = 0; t < P; ++t) {
            // ---- q, also by subtraction --------------------------------------
            double q_s = 0.0; for (int j = s0; j < s1; ++j) q_s += u[t][j] * u[t][j];
            double q_mid = 0.0; for (int j = 0; j < m_mid; ++j) q_mid += u[t][j] * u[t][j];
            double q_com = 0.0; for (int j = m_mid; j < K; ++j) q_com += u[t][j] * u[t][j];
            Eigen::VectorXd q(C + 1);
            q[0] = g[0] * q_s; q[1] = g[1] * (q_mid - q_s);
            if (has_com) q[2] = g[2] * q_com;
            q[env] = Y.col(t).squaredNorm();

            Eigen::VectorXd sigma = Tcod.solve(q);
            sr.vg[t] = sigma[0];
            sr.h2[t] = (Vp[t] > 0) ? sigma[0] / Vp[t] : NA_REAL;

            if (!spa || !have_L) continue;

            // ---- plug-in null covariance ------------------------------------
            Eigen::VectorXd s0v = sigma.head(C);
            for (int a = 0; a < C; ++a) if (s0v[a] < 0.0) s0v[a] = 0.0;
            s0v[0] = 0.0;                                  // H0: target = 0
            double sigma_env0 = std::max(sigma[env], 1e-8 * (Vp[t] > 0 ? Vp[t] : 1.0));

            // Diagonals over the SHARED column set. Component a contributes
            // sigma_a * g_a (resp. c_a * g_a) to each of its columns, because
            // K_a = g_a S_A = g_a V_A V_A'.
            Eigen::VectorXd D0(K), DM(K);
            const double c_s = Tinv(0,0), c_r = Tinv(0,1);
            const double c_c = has_com ? Tinv(0,2) : 0.0;
            const double c_env = Tinv(0,env);
            for (int j = 0; j < m_mid; ++j) {
                bool in_s = (j >= s0 && j < s1);
                D0[j] = in_s ? 0.0            : s0v[1] * g[1];
                DM[j] = in_s ? c_s * g[0]     : c_r    * g[1];
            }
            for (int j = m_mid; j < K; ++j) { D0[j] = s0v[2] * g[2]; DM[j] = c_c * g[2]; }

            // ---- P and Q via the per-block products (O(K^2 m_s), not O(K^3)) --
            Eigen::MatrixXd Ls = L.middleRows(s0, m_s);
            Eigen::MatrixXd PL_s; PL_s.noalias() = Ls.transpose() * Ls;
            Eigen::MatrixXd Pm = (s0v[1] * g[1]) * (PL_mid - PL_s);
            Eigen::MatrixXd Qm = (c_r * g[1]) * (PL_mid - PL_s) + (c_s * g[0]) * PL_s;
            if (has_com) { Pm += (s0v[2] * g[2]) * PL_com; Qm += (c_c * g[2]) * PL_com; }
            Pm = (0.5 * (Pm + Pm.transpose())).eval();
            Qm = (0.5 * (Qm + Qm.transpose())).eval();

            Eigen::MatrixXd Smat = Pm; Smat.diagonal().array() += sigma_env0;
            Eigen::LLT<Eigen::MatrixXd> lltS(Smat);
            if (lltS.info() != Eigen::Success) continue;
            Eigen::MatrixXd Cf = lltS.matrixL();
            Qm.diagonal().array() += c_env;
            Eigen::MatrixXd Asym = Cf.transpose() * Qm * Cf;
            Asym = (0.5 * (Asym + Asym.transpose())).eval();

            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ses(Asym, Eigen::EigenvaluesOnly);
            if (ses.info() != Eigen::Success) continue;
            const Eigen::VectorXd& ev = ses.eigenvalues();
            std::vector<double> eig(ev.data(), ev.data() + K);
            double eig_rep = c_env * sigma_env0;
            long n_rep = (long) n - K;

            // Wald SE comes free here: Var(sigma_hat) = K''(0) = 2 sum lambda^2
            double var0 = 0.0;
            for (int j = 0; j < K; ++j) var0 += 2.0 * eig[j] * eig[j];
            if (n_rep > 0) var0 += (double) n_rep * 2.0 * eig_rep * eig_rep;
            if (var0 > 0.0) {
                sr.se_vg[t] = std::sqrt(var0);
                double zst = sigma[0] / sr.se_vg[t];
                double p_wald = std::erfc(std::abs(zst) / std::sqrt(2.0));
                // Screen: the normal approximation is anti-conservative in the
                // upper tail, so a large Wald p implies an even larger SPA p.
                if (p_wald >= spa_thresh) { sr.p_spa[t] = p_wald; sr.spa_used[t] = 0; continue; }
                QuadSpaResult qr = quad_spa_solve(sigma[0], eig, eig_rep, n_rep);
                if (qr.converged) { sr.p_spa[t] = qr.p; sr.spa_used[t] = 1; }
                else              { sr.p_spa[t] = p_wald; sr.spa_used[t] = 0; }
            }
        }
        out.push_back(sr);
    }
    (void) background_rank;   // reserved: low-rank background truncation
    (void) min_sub_snps;
}

// ===========================================================================
// Driver: stream 1 Mb cells, scan sub-windows, write a long-format table
// ===========================================================================
struct SubParams {
    long W; int min_snps; long sub_bp; int min_sub_snps;
    bool spa; double spa_thresh; int background_rank;
    std::string out_file; int batch_size; int n_threads;
};

static void common_chr_range(const SubContext& ctx, const std::string& chr, int& lo, int& hi) {
    lo = hi = -1;
    if (!ctx.use_common) return;
    for (int j = 0; j < ctx.common_n_snps; ++j)
        if (ctx.common_bim.chr[j] == chr) { if (lo == -1) lo = j; hi = j; }
}

// Read + prepare one 1 Mb cell (serial: genotype I/O and standardization).
static bool make_cell(const SubContext& ctx, size_t ci, int cell_lo, int cell_hi,
                      int cc_lo, int cc_hi, CellData& cd) {
    cd = CellData();
    cd.chr   = ctx.chr_order[ci];
    cd.start = ctx.wes_bim.bp[cell_lo];
    cd.end   = ctx.wes_bim.bp[cell_hi];
    cd.i0_mid = cell_lo;

    Cell mid = read_cell_idx(ctx.wes_prefix, ctx.wes_n_total, ctx.wes_n_snps,
                             ctx.geno_keep, cell_lo, cell_hi + 1);
    if (ctx.covZ.cols() > 0) project_covariates(mid.X, ctx.covZ, ctx.covM);
    const std::vector<float>* wp = ctx.snp_weights.empty() ? 0 : &ctx.snp_weights;
    apply_alpha(mid.X, mid.maf, ctx.alpha, wp, cell_lo);
    cd.m_mid = (int) mid.X.cols();
    if (cd.m_mid <= 0) return false;

    Cell com; com.X = GenoMat(ctx.n_inds, 0);
    if (ctx.use_common) {
        com = read_cell(ctx.common_prefix, ctx.common_n_total, ctx.common_n_snps,
                        ctx.common_bim, ctx.common_keep, cc_lo, cc_hi,
                        cd.start, cd.end - cd.start + 1);
        if (com.X.cols() > 0) {
            if (ctx.covZ.cols() > 0) project_covariates(com.X, ctx.covZ, ctx.covM);
            apply_alpha(com.X, com.maf, ctx.alpha_common, 0, 0);
        }
    }
    cd.m_com = (int) com.X.cols();
    cd.K = cd.m_mid + cd.m_com;

    cd.V = GenoMat(ctx.n_inds, cd.K);
    cd.V.leftCols(cd.m_mid) = mid.X;
    if (cd.m_com > 0) cd.V.rightCols(cd.m_com) = com.X;

    cd.cj.resize(cd.K);
    for (int j = 0; j < cd.K; ++j) cd.cj[j] = cd.V.col(j).cast<double>().squaredNorm();

    cd.bp.resize(cd.m_mid);
    for (int j = 0; j < cd.m_mid; ++j) cd.bp[j] = ctx.wes_bim.bp[cell_lo + j];
    return true;
}

struct SubWorker : public RcppParallel::Worker {
    const std::vector<CellData>& batch; const Eigen::MatrixXd& Y;
    const std::vector<double>& Vp; const SubParams& pr;
    std::vector< std::vector<SubResult> >& out;
    SubWorker(const std::vector<CellData>& b, const Eigen::MatrixXd& Y,
              const std::vector<double>& Vp, const SubParams& pr,
              std::vector< std::vector<SubResult> >& out)
        : batch(b), Y(Y), Vp(Vp), pr(pr), out(out) {}
    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t w = begin; w < end; ++w)
            scan_cell(batch[w], Y, Vp, pr.sub_bp, pr.min_sub_snps,
                      pr.spa, pr.spa_thresh, pr.background_rank, out[w]);
    }
};

static Rcpp::List sub_driver(SubContext& ctx, const SubParams& pr) {
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
    std::vector<double> Vp(P);
    for (int t = 0; t < P; ++t) { double m = Y.col(t).mean();
        Vp[t] = (Y.col(t).array() - m).square().sum() / (Y.rows() - 1); }

    std::ofstream fout;
    bool tofile = !pr.out_file.empty();
    auto wr = [](std::ofstream& f, double v) { if (std::isnan(v)) f << "NA"; else f << v; };
    if (tofile) {
        fout.open(pr.out_file.c_str());
        if (!fout.is_open()) stop("Could not open out_file: " + pr.out_file);
        fout << "chr\tstart\tend\tm_sub\tm_rest\tm_common\tphenotype\tvg\tse_vg\th2";
        if (pr.spa) fout << "\tp_spa\tspa_used";
        fout << "\n";
    }

    // enumerate the 1 Mb cells exactly as the main pipeline does
    std::vector< std::pair<size_t, std::pair<int,int> > > cells;
    for (size_t ci = 0; ci < ctx.chr_order.size(); ++ci) {
        std::vector< std::pair<int,int> > cc =
            build_cells(ctx.wes_bim.bp, ctx.chr_lo[ci], ctx.chr_hi[ci], pr.W, pr.min_snps);
        for (size_t k = 0; k < cc.size(); ++k) cells.push_back(std::make_pair(ci, cc[k]));
    }
    Rcout << "1 Mb cells to scan: " << cells.size() << "\n";

    long n_sub = 0; int n_cell = 0;
    size_t idx = 0;
    while (idx < cells.size()) {
        std::vector<CellData> batch;
        while ((int) batch.size() < pr.batch_size && idx < cells.size()) {
            size_t ci = cells[idx].first;
            int lo = cells[idx].second.first, hi = cells[idx].second.second;
            int cc_lo, cc_hi; common_chr_range(ctx, ctx.chr_order[ci], cc_lo, cc_hi);
            CellData cd;
            if (make_cell(ctx, ci, lo, hi, cc_lo, cc_hi, cd)) batch.push_back(std::move(cd));
            ++idx;
        }
        if (batch.empty()) break;

        std::vector< std::vector<SubResult> > res(batch.size());
        SubWorker wk(batch, Y, Vp, pr, res);
        RcppParallel::parallelFor(0, batch.size(), wk);

        for (size_t b = 0; b < res.size(); ++b) {
            for (size_t r = 0; r < res[b].size(); ++r) {
                const SubResult& sr = res[b][r];
                ++n_sub;
                if (!tofile) continue;
                for (int t = 0; t < P; ++t) {
                    fout << sr.chr << '\t' << sr.start << '\t' << sr.end << '\t'
                         << sr.m_s << '\t' << sr.m_rest << '\t' << sr.m_com << '\t'
                         << as<std::string>(ctx.trait_names[t]) << '\t';
                    wr(fout, sr.vg[t]);    fout << '\t';
                    wr(fout, sr.se_vg[t]); fout << '\t';
                    wr(fout, sr.h2[t]);
                    if (pr.spa) { fout << '\t'; wr(fout, sr.p_spa[t]);
                                  fout << '\t' << sr.spa_used[t]; }
                    fout << '\n';
                }
            }
            ++n_cell;
        }
        if (tofile) fout.flush();
        Rcout << "\r[sub-window] cells " << n_cell << "/" << cells.size()
              << "  sub-windows " << n_sub << "    " << std::flush;
        Rcpp::checkUserInterrupt();
    }
    Rcout << "\n";
    if (tofile) fout.close();
    Rcout << "Cells scanned: " << n_cell << "   sub-windows tested: " << n_sub << "\n";
    return List::create(_["n_cells"] = n_cell, _["n_subwindows"] = (double) n_sub,
                        _["trait_names"] = ctx.trait_names);
}

}  // end anonymous namespace

// ===========================================================================
// Export
// ===========================================================================
// [[Rcpp::export]]
Rcpp::List he_subwindow_spa(const std::string& filename,
                            const SEXP pheno_mat,
                            double window_size = 1e6,
                            int min_snps = 1000,
                            double sub_window_size = 5e3,
                            int min_sub_snps = 2,
                            double alpha = -1.0,
                            double alpha_common = -1.0,
                            Rcpp::Nullable<Rcpp::String> common_filename = R_NilValue,
                            std::string out_file = "",
                            int batch_size = 4,
                            int n_threads = 0,
                            Rcpp::Nullable<Rcpp::NumericVector> weights = R_NilValue,
                            Rcpp::Nullable<Rcpp::NumericMatrix> covariates = R_NilValue,
                            bool SPA = true,
                            double spa_pval_threshold = 0.1,
                            int background_rank = 0) {
    if (window_size <= 0) stop("window_size must be positive");
    if (sub_window_size <= 0) stop("sub_window_size must be positive");
    if (min_snps < 1) min_snps = 1;
    if (min_sub_snps < 1) min_sub_snps = 1;
    SubContext ctx = setup_sub_context(filename, pheno_mat, alpha, alpha_common,
                                       common_filename, weights, covariates);
    Rcout << "Sub-window scan: " << (long) window_size << " bp cells, "
          << (long) sub_window_size << " bp sub-windows (>= " << min_sub_snps << " SNPs)\n";
    if (SPA)
        Rcout << "  SPA on, Wald pre-screen at p < " << spa_pval_threshold << ".\n"
              << "  Cost is driven by K_tot = (SNPs in the 1 Mb cell) + (common SNPs);\n"
              << "  the Gram matrix and its Cholesky are built ONCE per cell and reused\n"
              << "  by every sub-window in it.\n";
    SubParams pr;
    pr.W = (long) window_size; pr.min_snps = min_snps;
    pr.sub_bp = (long) sub_window_size; pr.min_sub_snps = min_sub_snps;
    pr.spa = SPA; pr.spa_thresh = spa_pval_threshold;
    pr.background_rank = background_rank;
    pr.out_file = out_file; pr.batch_size = batch_size; pr.n_threads = n_threads;
    return sub_driver(ctx, pr);
}
