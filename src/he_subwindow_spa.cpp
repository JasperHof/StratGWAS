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
    // Optional annotation: col 1 = flank flag, cols 2+ = middle categories.
    // snp_cat[j] = 0..n_annot_cat-1 (a tested category) or -1 (flank/background).
    bool use_annot; int n_annot_cat;
    std::vector<std::string> annot_names;
    std::vector<int> snp_cat;
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

    // Annotation: col 0 = flank flag, cols 1.. = middle functional categories.
    if (annotation.isNotNull()) {
        Rcpp::IntegerMatrix am(annotation.get());
        if (am.nrow() != ctx.wes_n_snps)
            stop("annotation must have one row per SNP in the WES .bim file");
        int ncol = am.ncol();
        if (ncol < 2) stop("annotation needs >= 2 columns: col 1 = flank flag, cols 2+ = middle categories");
        ctx.n_annot_cat = ncol - 1;
        ctx.snp_cat.assign(ctx.wes_n_snps, -1);
        for (int j = 0; j < ctx.wes_n_snps; ++j) {
            if (am(j, 0) == 1) { ctx.snp_cat[j] = -1; continue; }   // flank -> background
            for (int c = 1; c < ncol; ++c) if (am(j, c) == 1) { ctx.snp_cat[j] = c - 1; break; }
        }
        if (annot_names.isNotNull()) {
            Rcpp::CharacterVector nm(annot_names.get());
            if (nm.size() != ctx.n_annot_cat) stop("annot_names length must equal the number of category columns");
            for (int c = 0; c < ctx.n_annot_cat; ++c) ctx.annot_names.push_back(as<std::string>(nm[c]));
        } else {
            for (int c = 0; c < ctx.n_annot_cat; ++c) ctx.annot_names.push_back("cat" + std::to_string(c + 1));
        }
        ctx.use_annot = true;
        Rcout << "Annotation: " << ctx.n_annot_cat
              << " middle categories (+ flank column); an SPA p-value is reported per category\n";
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
static void test_chunk(const ChunkData& cd, const Eigen::MatrixXd& Y,
                       const std::vector<double>& Vp, bool spa, double spa_thresh,
                       ChunkResult& cr) {
    const int n = (int) Y.rows(), P = (int) Y.cols();
    const int K = cd.K, m_t = cd.m_t, m_f = cd.m_f, m_c = cd.m_c;
    const bool has_c = (m_c > 0);
    const int C = has_c ? 3 : 2, env = C;

    cr.chr = cd.chr; cr.start = cd.start; cr.end = cd.end;
    cr.m_t = m_t; cr.m_f = m_f; cr.m_c = m_c;
    cr.vg.assign(P, NA_REAL); cr.se_vg.assign(P, NA_REAL);
    cr.h2.assign(P, NA_REAL); cr.p_spa.assign(P, NA_REAL);
    cr.spa_used.assign(P, 0);
    if (m_t <= 0 || m_f <= 0) return;

    // ---- Gram matrix (K x K, K is BOUNDED by chunk_size and flank_chunks) ---
    GenoMat Gf = GenoMat::Zero(K, K);
    Gf.selfadjointView<Eigen::Upper>().rankUpdate(cd.V.transpose());
    Eigen::MatrixXd G = GenoMat(Gf.selfadjointView<Eigen::Upper>()).cast<double>();
    Eigen::MatrixXd G2 = G.array().square();

    const int o_t = 0, o_f = m_t, o_c = m_t + m_f;
    double tr_t = 0, tr_f = 0, tr_c = 0;
    for (int j = 0; j < m_t; ++j) tr_t += cd.cj[o_t + j];
    for (int j = 0; j < m_f; ++j) tr_f += cd.cj[o_f + j];
    for (int j = 0; j < m_c; ++j) tr_c += cd.cj[o_c + j];
    if (!(tr_t > 0.0) || !(tr_f > 0.0)) return;

    std::vector<double> g(C);
    g[0] = (double) n / tr_t; g[1] = (double) n / tr_f;
    if (has_c) { if (!(tr_c > 0.0)) return; g[2] = (double) n / tr_c; }

    // ---- moment matrix T:  tr(S_a S_b) = ||G[A,B]||_F^2 ---------------------
    Eigen::MatrixXd T = Eigen::MatrixXd::Zero(C + 1, C + 1);
    T(0,0) = g[0]*g[0]*G2.block(o_t,o_t,m_t,m_t).sum();
    T(1,1) = g[1]*g[1]*G2.block(o_f,o_f,m_f,m_f).sum();
    T(0,1) = T(1,0) = g[0]*g[1]*G2.block(o_t,o_f,m_t,m_f).sum();
    if (has_c) {
        T(2,2) = g[2]*g[2]*G2.block(o_c,o_c,m_c,m_c).sum();
        T(0,2) = T(2,0) = g[0]*g[2]*G2.block(o_t,o_c,m_t,m_c).sum();
        T(1,2) = T(2,1) = g[1]*g[2]*G2.block(o_f,o_c,m_f,m_c).sum();
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

    for (int t = 0; t < P; ++t) {
        Eigen::VectorXd u = (cd.V.transpose() * Y.col(t).cast<float>()).cast<double>();
        Eigen::VectorXd q(C + 1);
        double q_t = 0, q_f = 0, q_c = 0;
        for (int j = 0; j < m_t; ++j) q_t += u[o_t+j]*u[o_t+j];
        for (int j = 0; j < m_f; ++j) q_f += u[o_f+j]*u[o_f+j];
        for (int j = 0; j < m_c; ++j) q_c += u[o_c+j]*u[o_c+j];
        q[0] = g[0]*q_t; q[1] = g[1]*q_f; if (has_c) q[2] = g[2]*q_c;
        q[env] = Y.col(t).squaredNorm();

        Eigen::VectorXd sigma = Tcod.solve(q);
        cr.vg[t] = sigma[0];
        cr.h2[t] = (Vp[t] > 0) ? sigma[0] / Vp[t] : NA_REAL;
        if (!spa || !have_L) continue;

        // ---- plug-in null covariance ---------------------------------------
        Eigen::VectorXd s0v = sigma.head(C);
        for (int a = 0; a < C; ++a) if (s0v[a] < 0.0) s0v[a] = 0.0;
        s0v[0] = 0.0;
        double sigma_env0 = std::max(sigma[env], 1e-8 * (Vp[t] > 0 ? Vp[t] : 1.0));

        const double c_t = Tinv(0,0), c_f = Tinv(0,1);
        const double c_c = has_c ? Tinv(0,2) : 0.0, c_env = Tinv(0,env);
        Eigen::VectorXd D0(K), DM(K);
        for (int j = 0; j < m_t; ++j) { D0[o_t+j] = 0.0;              DM[o_t+j] = c_t * g[0]; }
        for (int j = 0; j < m_f; ++j) { D0[o_f+j] = s0v[1] * g[1];    DM[o_f+j] = c_f * g[1]; }
        for (int j = 0; j < m_c; ++j) { D0[o_c+j] = s0v[2] * g[2];    DM[o_c+j] = c_c * g[2]; }

        Eigen::MatrixXd Pm = L.transpose() * (D0.asDiagonal() * L);
        Eigen::MatrixXd Qm = L.transpose() * (DM.asDiagonal() * L);
        Pm = (0.5*(Pm + Pm.transpose())).eval();
        Qm = (0.5*(Qm + Qm.transpose())).eval();
        Eigen::MatrixXd Smat = Pm; Smat.diagonal().array() += sigma_env0;
        Eigen::LLT<Eigen::MatrixXd> lltS(Smat);
        if (lltS.info() != Eigen::Success) continue;
        Eigen::MatrixXd Cf = lltS.matrixL();
        Qm.diagonal().array() += c_env;
        Eigen::MatrixXd Asym = Cf.transpose() * Qm * Cf;
        Asym = (0.5*(Asym + Asym.transpose())).eval();

        // ---- WALD SCREEN, BEFORE THE EIGENSOLVE -----------------------------
        // The explicit eigenvalues are those of the symmetric Asym, so
        //     Var(sigma_hat) = 2 sum lambda^2 = 2 ||Asym||_F^2
        // exactly, in O(K^2). Computing it here (rather than after the
        // eigensolve, as an earlier version mistakenly did) is what makes the
        // screen worth anything: it skips the ~10 K^3 eigensolve entirely for
        // the ~90% of chunks that are nowhere near significant.
        double eig_rep = c_env * sigma_env0;
        long n_rep = (long) n - K;
        double var0 = 2.0 * Asym.squaredNorm();
        if (n_rep > 0) var0 += (double) n_rep * 2.0 * eig_rep * eig_rep;
        if (!(var0 > 0.0)) continue;
        cr.se_vg[t] = std::sqrt(var0);
        double p_wald = std::erfc(std::abs(sigma[0] / cr.se_vg[t]) / std::sqrt(2.0));
        if (p_wald >= spa_thresh) { cr.p_spa[t] = p_wald; cr.spa_used[t] = 0; continue; }

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
    long common_bp; int max_common_snps;
    bool spa; double spa_thresh;
    std::string out_file; int batch_size, n_threads;
};

static void common_chr_range(const ChunkContext& ctx, const std::string& chr, int& lo, int& hi) {
    lo = hi = -1;
    if (!ctx.use_common) return;
    for (int j = 0; j < ctx.common_n_snps; ++j)
        if (ctx.common_bim.chr[j] == chr) { if (lo == -1) lo = j; hi = j; }
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
                       int flank_snps, long common_bp, int max_common, ChunkData& cd) {
    cd = ChunkData();
    cd.chr = ctx.chr_order[ci];
    cd.start = ctx.wes_bim.bp[a]; cd.end = ctx.wes_bim.bp[b - 1];
    cd.i0_target = a;

    const int fL0 = std::max(chr_lo, a - flank_snps), fL1 = a;
    const int fR0 = b, fR1 = std::min(chr_hi + 1, b + flank_snps);

    const std::vector<float>* wp = ctx.snp_weights.empty() ? 0 : &ctx.snp_weights;
    Cell tgt = read_cell_idx(ctx.wes_prefix, ctx.wes_n_total, ctx.wes_n_snps, ctx.geno_keep, a, b);
    if (ctx.covZ.cols() > 0) project_covariates(tgt.X, ctx.covZ, ctx.covM);
    apply_alpha(tgt.X, tgt.maf, ctx.alpha, wp, a);
    cd.m_t = (int) tgt.X.cols();
    if (cd.m_t <= 0) return false;

    Cell fl; fl.X = GenoMat(ctx.n_inds, 0);
    Cell fr; fr.X = GenoMat(ctx.n_inds, 0);
    if (fL1 > fL0) { fl = read_cell_idx(ctx.wes_prefix, ctx.wes_n_total, ctx.wes_n_snps, ctx.geno_keep, fL0, fL1);
                     if (ctx.covZ.cols() > 0) project_covariates(fl.X, ctx.covZ, ctx.covM);
                     apply_alpha(fl.X, fl.maf, ctx.alpha, wp, fL0); }
    if (fR1 > fR0) { fr = read_cell_idx(ctx.wes_prefix, ctx.wes_n_total, ctx.wes_n_snps, ctx.geno_keep, fR0, fR1);
                     if (ctx.covZ.cols() > 0) project_covariates(fr.X, ctx.covZ, ctx.covM);
                     apply_alpha(fr.X, fr.maf, ctx.alpha, wp, fR0); }
    cd.m_f = (int) fl.X.cols() + (int) fr.X.cols();
    if (cd.m_f <= 0) return false;

    Cell com; com.X = GenoMat(ctx.n_inds, 0);
    if (ctx.use_common) {
        // FIXED bp window centred on the TARGET chunk -- NOT the span of the
        // chunk+flanks. SNP spacing in exome data is wildly uneven (consecutive
        // SNPs in this dataset can be ~800 kb apart), so the span of a fixed
        // NUMBER of SNPs is unbounded: a region straddling a couple of
        // intergenic gaps would drag in thousands of common SNPs and push K_tot
        // straight back to where the bp-window designs failed. A fixed window
        // plus a hard cap keeps K_tot predictable, which is the entire point of
        // the chunk design.
        long ctr = (cd.start + cd.end) / 2;
        long half = common_bp / 2;
        long lo_bp = (ctr > half) ? (ctr - half) : 0;
        com = read_cell(ctx.common_prefix, ctx.common_n_total, ctx.common_n_snps,
                        ctx.common_bim, ctx.common_keep, cc_lo, cc_hi, lo_bp, common_bp);
        thin_cell(com, max_common);
        if (com.X.cols() > 0) {
            if (ctx.covZ.cols() > 0) project_covariates(com.X, ctx.covZ, ctx.covM);
            apply_alpha(com.X, com.maf, ctx.alpha_common, 0, 0);
        }
    }
    cd.m_c = (int) com.X.cols();
    cd.K = cd.m_t + cd.m_f + cd.m_c;

    cd.V = GenoMat(ctx.n_inds, cd.K);
    cd.V.leftCols(cd.m_t) = tgt.X;
    int off = cd.m_t;
    if (fl.X.cols() > 0) { cd.V.middleCols(off, fl.X.cols()) = fl.X; off += (int) fl.X.cols(); }
    if (fr.X.cols() > 0) { cd.V.middleCols(off, fr.X.cols()) = fr.X; off += (int) fr.X.cols(); }
    if (cd.m_c > 0) cd.V.rightCols(cd.m_c) = com.X;

    cd.cj.resize(cd.K);
    for (int j = 0; j < cd.K; ++j) cd.cj[j] = cd.V.col(j).cast<double>().squaredNorm();
    return true;
}

struct ChunkWorker : public RcppParallel::Worker {
    const std::vector<ChunkData>& batch; const Eigen::MatrixXd& Y;
    const std::vector<double>& Vp; const ChunkParams& pr;
    std::vector<ChunkResult>& out;
    ChunkWorker(const std::vector<ChunkData>& b, const Eigen::MatrixXd& Y,
                const std::vector<double>& Vp, const ChunkParams& pr, std::vector<ChunkResult>& out)
        : batch(b), Y(Y), Vp(Vp), pr(pr), out(out) {}
    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t w = begin; w < end; ++w)
            test_chunk(batch[w], Y, Vp, pr.spa, pr.spa_thresh, out[w]);
    }
};

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
    std::vector<double> Vp(P);
    for (int t = 0; t < P; ++t) { double m = Y.col(t).mean();
        Vp[t] = (Y.col(t).array() - m).square().sum() / (Y.rows() - 1); }

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
    struct Job { size_t ci; int a, b, chr_lo, chr_hi; };
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
    Rcout << "  K_tot per chunk <= " << (pr.chunk_size * (1 + 2 * pr.flank_chunks) + pr.max_common_snps)
          << "  (" << (pr.chunk_size * (1 + 2 * pr.flank_chunks)) << " WES + <= "
          << pr.max_common_snps << " common)\n";
    if (ctx.use_common)
        Rcout << "  Common GRM: fixed " << pr.common_bp
              << " bp window centred on each chunk, thinned to <= "
              << pr.max_common_snps << " SNPs -- so K_tot is BOUNDED regardless\n"
              << "  of how far the chunk's SNPs happen to span.\n";

    long n_done = 0; size_t idx = 0;
    typedef std::chrono::steady_clock clk;
    clk::time_point t0 = clk::now();
    while (idx < jobs.size()) {
        std::vector<ChunkData> batch;
        while ((int) batch.size() < pr.batch_size && idx < jobs.size()) {
            const Job& j = jobs[idx];
            int cc_lo, cc_hi; common_chr_range(ctx, ctx.chr_order[j.ci], cc_lo, cc_hi);
            ChunkData cd;
            if (make_chunk(ctx, j.ci, j.a, j.b, j.chr_lo, j.chr_hi, cc_lo, cc_hi,
                           flank_snps, pr.common_bp, pr.max_common_snps, cd))
                batch.push_back(std::move(cd));
            ++idx;
        }
        if (batch.empty()) break;
        std::vector<ChunkResult> res(batch.size());
        ChunkWorker wk(batch, Y, Vp, pr, res);
        RcppParallel::parallelFor(0, batch.size(), wk);

        for (size_t b = 0; b < res.size(); ++b) {
            const ChunkResult& cr = res[b];
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
        double eta = (n_done > 0) ? el * ((double) jobs.size() - n_done) / n_done : 0.0;
        Rcout << "\r[chunk] " << n_done << "/" << jobs.size()
              << "  " << (int)(100.0 * n_done / jobs.size()) << "%"
              << "  elapsed " << (int) el << "s  eta " << (int) eta << "s     " << std::flush;
        Rcpp::checkUserInterrupt();
    }
    Rcout << "\n";
    if (tofile) fout.close();
    Rcout << "Chunks tested: " << n_done << "\n";
    return List::create(_["n_chunks"] = (double) n_done, _["trait_names"] = ctx.trait_names);
}

// ===========================================================================
// ANNOTATION PATH (only used when an annotation matrix is supplied).
//
// Instead of one target = the whole chunk, the chunk's middle SNPs are split
// by the annotation into functional CATEGORIES (annotation cols 2+). Column 1
// of the annotation flags flank/background SNPs. The model per chunk is:
//   { cat_0, cat_1, ..., cat_{A-1}  (each SPA-tested),
//     flank (neighbour chunks + any middle SNP flagged flank or uncategorized),
//     [common], sigma_e I }
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
                             int flank_snps, long common_bp, int max_common, ChunkDataA& cd) {
    cd = ChunkDataA();
    cd.chr = ctx.chr_order[ci];
    cd.start = ctx.wes_bim.bp[a]; cd.end = ctx.wes_bim.bp[b - 1];
    const std::vector<float>* wp = ctx.snp_weights.empty() ? 0 : &ctx.snp_weights;

    Cell tgt = read_cell_idx(ctx.wes_prefix, ctx.wes_n_total, ctx.wes_n_snps, ctx.geno_keep, a, b);
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
    if (fL1 > fL0) { fl = read_cell_idx(ctx.wes_prefix, ctx.wes_n_total, ctx.wes_n_snps, ctx.geno_keep, fL0, fL1);
                     if (ctx.covZ.cols() > 0) project_covariates(fl.X, ctx.covZ, ctx.covM);
                     apply_alpha(fl.X, fl.maf, ctx.alpha, wp, fL0); }
    if (fR1 > fR0) { fr = read_cell_idx(ctx.wes_prefix, ctx.wes_n_total, ctx.wes_n_snps, ctx.geno_keep, fR0, fR1);
                     if (ctx.covZ.cols() > 0) project_covariates(fr.X, ctx.covZ, ctx.covM);
                     apply_alpha(fr.X, fr.maf, ctx.alpha, wp, fR0); }
    int m_flank = (int) bg_cols.size() + (int) fl.X.cols() + (int) fr.X.cols();
    if (m_flank <= 0) return false;                // model needs a background component

    // common GRM (fixed bp window centred on the chunk, thinned) -- as no-annot
    Cell com; com.X = GenoMat(ctx.n_inds, 0);
    if (ctx.use_common) {
        long ctr = (cd.start + cd.end) / 2, half = common_bp / 2;
        long lo_bp = (ctr > half) ? (ctr - half) : 0;
        com = read_cell(ctx.common_prefix, ctx.common_n_total, ctx.common_n_snps,
                        ctx.common_bim, ctx.common_keep, cc_lo, cc_hi, lo_bp, common_bp);
        thin_cell(com, max_common);
        if (com.X.cols() > 0) {
            if (ctx.covZ.cols() > 0) project_covariates(com.X, ctx.covZ, ctx.covM);
            apply_alpha(com.X, com.maf, ctx.alpha_common, 0, 0);
        }
    }
    int m_c = (int) com.X.cols();

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
    for (size_t q = 0; q < bg_cols.size(); ++q) cd.V.col(off++) = tgt.X.col(bg_cols[q]);
    if (fl.X.cols() > 0) { cd.V.middleCols(off, fl.X.cols()) = fl.X; off += (int) fl.X.cols(); }
    if (fr.X.cols() > 0) { cd.V.middleCols(off, fr.X.cols()) = fr.X; off += (int) fr.X.cols(); }
    if (m_c > 0) cd.V.rightCols(m_c) = com.X;
    cd.m_flank = m_flank; cd.m_common = m_c;
    cd.cj.resize(cd.K);
    for (int j = 0; j < cd.K; ++j) cd.cj[j] = cd.V.col(j).cast<double>().squaredNorm();
    return true;
}

// Multi-category SPA test. Components: 0..A-1 categories (tested), A flank,
// A+1 common (if any), env residual. Mirrors test_chunk but loops the SPA over
// the A categories, reusing the single per-chunk Cholesky.
static void test_chunk_annot(const ChunkDataA& cd, const Eigen::MatrixXd& Y,
                             const std::vector<double>& Vp, bool spa, double spa_thresh,
                             ChunkResultA& cr) {
    const int n = (int) Y.rows(), P = (int) Y.cols(), K = cd.K;
    const int A = (int) cd.cat_m.size();
    const bool has_c = (cd.m_common > 0);
    const int C = A + 1 + (has_c ? 1 : 0), env = C;   // +1 flank, +1 common

    cr.chr = cd.chr; cr.start = cd.start; cr.end = cd.end;
    cr.m_flank = cd.m_flank; cr.m_common = cd.m_common;
    cr.cat_name = cd.cat_name; cr.cat_m = cd.cat_m;
    cr.vg.assign(A, std::vector<double>(P, NA_REAL));
    cr.se_vg.assign(A, std::vector<double>(P, NA_REAL));
    cr.h2.assign(A, std::vector<double>(P, NA_REAL));
    cr.p_spa.assign(A, std::vector<double>(P, NA_REAL));
    cr.spa_used.assign(A, std::vector<int>(P, 0));
    if (A <= 0 || cd.m_flank <= 0) return;

    // component column offsets in V: cats, then flank, then common
    std::vector<int> off(C + 1, 0);
    for (int c = 0; c < A; ++c) off[c + 1] = off[c] + cd.cat_m[c];
    off[A + 1] = off[A] + cd.m_flank;
    if (has_c) off[A + 2] = off[A + 1] + cd.m_common;

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

    for (int t = 0; t < P; ++t) {
        Eigen::VectorXd u = (cd.V.transpose() * Y.col(t).cast<float>()).cast<double>();
        Eigen::VectorXd q(C + 1);
        for (int c = 0; c < C; ++c) {
            double qc = 0.0; for (int j = off[c]; j < off[c + 1]; ++j) qc += u[j] * u[j];
            q[c] = g[c] * qc;
        }
        q[env] = Y.col(t).squaredNorm();
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
            Eigen::VectorXd s0v = s0base; s0v[c] = 0.0;   // null the tested category
            Eigen::VectorXd D0(K), DM(K);
            const double c_env = Tinv(c, env);
            for (int a = 0; a < C; ++a) {
                double d0 = (a == c) ? 0.0 : s0v[a] * g[a];
                double dm = Tinv(c, a) * g[a];
                for (int j = off[a]; j < off[a + 1]; ++j) { D0[j] = d0; DM[j] = dm; }
            }
            Eigen::MatrixXd Pm = L.transpose() * (D0.asDiagonal() * L);
            Eigen::MatrixXd Qm = L.transpose() * (DM.asDiagonal() * L);
            Pm = (0.5 * (Pm + Pm.transpose())).eval();
            Qm = (0.5 * (Qm + Qm.transpose())).eval();
            Eigen::MatrixXd Smat = Pm; Smat.diagonal().array() += sigma_env0;
            Eigen::LLT<Eigen::MatrixXd> lltS(Smat);
            if (lltS.info() != Eigen::Success) continue;
            Eigen::MatrixXd Cf = lltS.matrixL();
            Qm.diagonal().array() += c_env;
            Eigen::MatrixXd Asym = Cf.transpose() * Qm * Cf;
            Asym = (0.5 * (Asym + Asym.transpose())).eval();

            double eig_rep = c_env * sigma_env0;
            long n_rep = (long) n - K;
            double var0 = 2.0 * Asym.squaredNorm();
            if (n_rep > 0) var0 += (double) n_rep * 2.0 * eig_rep * eig_rep;
            if (!(var0 > 0.0)) continue;
            cr.se_vg[c][t] = std::sqrt(var0);
            double p_wald = std::erfc(std::abs(sigma[c] / cr.se_vg[c][t]) / std::sqrt(2.0));
            if (p_wald >= spa_thresh) { cr.p_spa[c][t] = p_wald; cr.spa_used[c][t] = 0; continue; }

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
    const std::vector<ChunkDataA>& batch; const Eigen::MatrixXd& Y;
    const std::vector<double>& Vp; const ChunkParams& pr;
    std::vector<ChunkResultA>& out;
    ChunkWorkerA(const std::vector<ChunkDataA>& b, const Eigen::MatrixXd& Y,
                 const std::vector<double>& Vp, const ChunkParams& pr, std::vector<ChunkResultA>& out)
        : batch(b), Y(Y), Vp(Vp), pr(pr), out(out) {}
    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t w = begin; w < end; ++w)
            test_chunk_annot(batch[w], Y, Vp, pr.spa, pr.spa_thresh, out[w]);
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
    std::vector<double> Vp(P);
    for (int t = 0; t < P; ++t) { double m = Y.col(t).mean();
        Vp[t] = (Y.col(t).array() - m).square().sum() / (Y.rows() - 1); }

    std::ofstream fout; bool tofile = !pr.out_file.empty();
    auto wr = [](std::ofstream& f, double v) { if (std::isnan(v)) f << "NA"; else f << v; };
    if (tofile) {
        fout.open(pr.out_file.c_str());
        if (!fout.is_open()) stop("Could not open out_file: " + pr.out_file);
        fout << "chr\tstart\tend\tcategory\tm_cat\tm_flank\tm_common\tphenotype\tvg\tse_vg\th2";
        if (pr.spa) fout << "\tp_spa\tspa_used";
        fout << "\n";
    }

    struct Job { size_t ci; int a, b, chr_lo, chr_hi; };
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

    long n_done = 0; size_t idx = 0;
    typedef std::chrono::steady_clock clk;
    clk::time_point t0 = clk::now();
    while (idx < jobs.size()) {
        std::vector<ChunkDataA> batch;
        while ((int) batch.size() < pr.batch_size && idx < jobs.size()) {
            const Job& j = jobs[idx];
            int cc_lo, cc_hi; common_chr_range(ctx, ctx.chr_order[j.ci], cc_lo, cc_hi);
            ChunkDataA cd;
            if (make_chunk_annot(ctx, j.ci, j.a, j.b, j.chr_lo, j.chr_hi, cc_lo, cc_hi,
                                 flank_snps, pr.common_bp, pr.max_common_snps, cd))
                batch.push_back(std::move(cd));
            ++idx;
        }
        if (batch.empty()) break;
        std::vector<ChunkResultA> res(batch.size());
        ChunkWorkerA wk(batch, Y, Vp, pr, res);
        RcppParallel::parallelFor(0, batch.size(), wk);

        for (size_t bb = 0; bb < res.size(); ++bb) {
            const ChunkResultA& cr = res[bb];
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
        double eta = (n_done > 0) ? el * ((double) jobs.size() - n_done) / n_done : 0.0;
        Rcout << "\r[chunk-annot] " << n_done << "/" << jobs.size()
              << "  " << (int)(100.0 * n_done / jobs.size()) << "%"
              << "  elapsed " << (int) el << "s  eta " << (int) eta << "s     " << std::flush;
        Rcpp::checkUserInterrupt();
    }
    Rcout << "\n";
    if (tofile) fout.close();
    Rcout << "Chunks tested: " << n_done << "\n";
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
                        double common_window = 1e6,
                        int max_common_snps = 400,
                        std::string out_file = "",
                        int batch_size = 16,
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
    pr.common_bp = (long) common_window; pr.max_common_snps = max_common_snps;
    pr.spa = SPA; pr.spa_thresh = spa_pval_threshold;
    pr.out_file = out_file; pr.batch_size = batch_size; pr.n_threads = n_threads;
    // Annotation supplied -> per-category SPA path. There the flank/background
    // comes from the annotation (column 1 + uncategorized middle SNPs), so
    // flank_chunks = 0 is allowed (no neighbouring-chunk flank is read).
    if (ctx.use_annot) return chunk_driver_annot(ctx, pr);
    // No annotation: the neighbouring chunks are the ONLY background, so at
    // least one flank chunk is required.
    if (flank_chunks < 1) stop("flank_chunks must be >= 1 when no annotation is given (the model needs a background component)");
    return chunk_driver(ctx, pr);
}
