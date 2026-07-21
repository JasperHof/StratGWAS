// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppParallel)]]
#include <RcppEigen.h>
#include <Rcpp.h>
#include <RcppParallel.h>
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
// FLEX-style PARTITIONED sliding-window heritability.
//
// Same sliding-window / flanking design as he_sliding_window.cpp, but the
// MIDDLE window is split into user-supplied functional CATEGORIES (e.g. MAF
// bins), each fit as its own variance component. Left/right WES flanks and an
// optional common-SNP GRM remain single background components. We record the
// per-category h2 of the middle window and slide.
//
//   components per window = { left, right, [common] }  (background)
//                         + { cat_1, cat_2, ... }      (middle, RECORDED)
//                         + residual (I)
//
// WINDOW SIZING (window_size + min_snps): the genome is partitioned into
// consecutive CELLS, each of which is at least `window_size` bp AND at least
// `min_snps` SNPs. A cell grows past window_size until it has min_snps SNPs (so
// sparse regions are not chopped into tiny, high-variance windows). Each cell
// serves as a middle exactly once and as a flank for its two neighbours, so the
// rolling cache still reads every cell only once (right -> middle -> left).
//
// Two knobs beyond he_sliding_window:
//   * snp_cat / cat_names : an (n_snps x n_cat) 0/1 membership matrix aligned to
//     the WES .bim, defining the categories. Only MIDDLE-window SNPs are
//     partitioned; middle SNPs in no category go into an "uncategorized" bin.
//   * alpha (default -1): each standardized column is reweighted by
//         w = [2 f (1-f)]^{(1+alpha)/2}      (f = MAF),
//     the LDAK/GRE heritability-model knob. alpha = -1 -> w = 1 (unweighted).
//
// IMPORTANT (why alpha needs trace normalization): the variance component s_c
// equals the heritability contribution ONLY when tr(G_c) = n (mean GRM diagonal
// 1). alpha-weighting changes tr(G_c), so we scale each component's columns so
// that tr(G_c) = n; then s_c is directly the (per-category) h2 in phenotype-
// variance units, for ANY alpha. (Skipping this makes estimates correct only at
// alpha = -1.)
//
// Estimators (same as he_sliding_window):
//   he_sliding_window_part()   - randomized Haseman-Elston (method of moments).
//   reml_sliding_window_part()  - low-rank / Woodbury AI-REML.
// Both parallelize over windows and stream a long-format per-category table.
// ===========================================================================

// ---------------------------------------------------------------------------
// EVERYTHING below lives in an anonymous namespace so its types/functions have
// INTERNAL linkage and cannot collide with identically-named symbols in the
// other .cpp files of the package (e.g. `struct RunContext`, `setup_context`,
// which he_sliding_window.cpp also defines but with a DIFFERENT layout). Two
// differing definitions of the same externally-visible name across translation
// units is an ODR violation -> undefined behaviour -> the segfault you get only
// when the whole package is compiled together (sourceCpp compiles this file
// alone, so it never triggers). The two [[Rcpp::export]] wrappers sit OUTSIDE
// this namespace (they need external linkage) and call into it.
// ---------------------------------------------------------------------------
namespace {

typedef Eigen::MatrixXf GenoMat;

// --------------------------------------------------------------------------
// Shared low-level helpers (same conventions as he_sliding_window.cpp)
// --------------------------------------------------------------------------
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
    std::vector<float> maf;        // per-column folded MAF
    int i0;                        // .bim index of column 0 (-1 if empty)
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
                          const std::vector<int>& keep, int i0, int i1) {
    Cell cell; cell.i0 = -1;
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

// Build ONE variance-component matrix from a set of source columns: apply the
// alpha weight per column (times the optional per-SNP LD weight), then scale the
// whole block so tr(G) = n (so the fitted variance component is directly the h2
// contribution). Returns an empty (0-col) matrix if the block is degenerate.
//
// LD WEIGHTS: `wts` (may be NULL) is the per-SNP weight vector aligned to the WES
// .bim, and `i0` is the .bim index of column 0 of X, so the weight of local column
// j is wts[i0 + j]. A SNP's variance contribution scales with its weight, so the
// COLUMN multiplier is sqrt(weight). The trace normalization below is unchanged,
// so the fitted component is still directly the h2 contribution.
static GenoMat build_component(const GenoMat& X, const std::vector<float>& maf,
                               const std::vector<int>& cols, double alpha, int n,
                               const std::vector<float>* wts = 0, int i0 = 0) {
    int m = (int) cols.size();
    if (m == 0) return GenoMat(X.rows(), 0);
    GenoMat W(X.rows(), m);
    bool unit = (alpha == -1.0);                 // alpha weight = 1 everywhere; skip pow
    double e = (1.0 + alpha) / 2.0;
    int nw = (wts != 0) ? (int) wts->size() : 0;
    for (int k = 0; k < m; ++k) {
        int j = cols[k];
        double w = 1.0;
        if (!unit) {
            double f = maf[j];
            w = (f > 0.0 && f < 1.0) ? std::pow(2.0 * f * (1.0 - f), e) : 1.0;
        }
        if (nw > 0) {                            // per-SNP LD weight (sqrt on the column)
            int g = i0 + j;
            double lw = (g >= 0 && g < nw) ? (double) (*wts)[g] : 1.0;
            if (lw < 0.0) lw = 0.0;
            w *= std::sqrt(lw);
        }
        if (w == 1.0) { W.col(k) = X.col(j); continue; }
        W.col(k) = X.col(j) * (float) w;
    }
    double ss = W.colwise().squaredNorm().cast<double>().sum();   // ||W||_F^2, double
    if (ss < 1e-12) return GenoMat(X.rows(), 0);                  // all-zero -> drop
    double c = std::sqrt((double) n * m / ss);                   // tr(G)=n after scaling
    W *= (float) c;
    return W;
}

// --------------------------------------------------------------------------
// Context
// --------------------------------------------------------------------------
struct RunContext {
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
    // categories
    int n_cat;
    std::vector<std::string> cat_names;      // length n_cat (+ "uncategorized" appended)
    std::vector< std::vector<int> > snp_cats; // per WES SNP: category ids it belongs to
    double alpha;
    // Optional per-SNP LD weights, aligned to the WES .bim (empty = all weights 1).
    // Applied to the WES components (flanks + middle categories) only; the
    // common-SNP background is a different fileset so it stays unweighted.
    std::vector<float> snp_weights;
};

static RunContext setup_context(const std::string& filename, const SEXP pheno_mat,
                                const IntegerMatrix& snp_cat, const CharacterVector& cat_names,
                                double alpha, Rcpp::Nullable<Rcpp::String> common_filename,
                                Rcpp::Nullable<Rcpp::NumericVector> weights) {
    RunContext ctx;
    ctx.wes_prefix = filename;
    ctx.wes_n_snps = count_lines(filename + ".bim");
    List fam = read_fam_file(filename);
    CharacterVector geno_iid = fam["iid"];
    ctx.wes_n_total = geno_iid.size();
    ctx.alpha = alpha;
    Rcout << "WES individuals (.fam): " << ctx.wes_n_total << "\n";
    Rcout << "WES SNPs (.bim): "        << ctx.wes_n_snps  << "\n";

    ctx.wes_bim = read_bim_positions(filename + ".bim");
    if (ctx.wes_bim.n_snps != ctx.wes_n_snps)
        stop("Mismatch between .bim line count and parsed positions");

    if (snp_cat.nrow() != ctx.wes_n_snps)
        stop("snp_cat must have one row per SNP in the WES .bim file");
    if (snp_cat.ncol() != cat_names.size())
        stop("snp_cat must have one column per entry in cat_names");
    ctx.n_cat = cat_names.size();
    for (int c = 0; c < ctx.n_cat; ++c) ctx.cat_names.push_back(as<std::string>(cat_names[c]));
    ctx.cat_names.push_back("uncategorized");     // catch-all = category index n_cat
    // per-SNP membership (thread-safe plain C++)
    ctx.snp_cats.resize(ctx.wes_n_snps);
    for (int j = 0; j < ctx.wes_n_snps; ++j)
        for (int c = 0; c < ctx.n_cat; ++c)
            if (snp_cat(j, c) == 1) ctx.snp_cats[j].push_back(c);

    // Optional per-SNP LD weights (one per WES .bim row). Stored as plain floats
    // so the parallel workers can read them without touching the R API.
    if (weights.isNotNull()) {
        Rcpp::NumericVector wv(weights.get());
        if (wv.size() != ctx.wes_n_snps)
            stop("weights must have one entry per SNP in the WES .bim file");
        ctx.snp_weights.resize(ctx.wes_n_snps);
        for (int j = 0; j < ctx.wes_n_snps; ++j) {
            double w = wv[j];
            if (ISNAN(w) || w < 0.0) w = 0.0;      // NA / negative -> drop the SNP
            ctx.snp_weights[j] = (float) w;
        }
        Rcout << "Using per-SNP LD weights (" << ctx.wes_n_snps << " values)\n";
    }

    // Phenotype
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
        for (int j = 0; j < ctx.n_pheno; ++j)
            ctx.Y(i, j) = pheno(pheno_keep[i], j);

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

static void common_chr_range(const RunContext& ctx, const std::string& chr, int& cc_lo, int& cc_hi) {
    cc_lo = cc_hi = -1;
    if (!ctx.use_common) return;
    for (int j = 0; j < ctx.common_n_snps; ++j)
        if (ctx.common_bim.chr[j] == chr) { if (cc_lo == -1) cc_lo = j; cc_hi = j; }
}

// --------------------------------------------------------------------------
// One assembled window: all component matrices (background + category), with
// bookkeeping of which components are the recorded (middle-category) signal.
// --------------------------------------------------------------------------
struct PartWindow {
    std::string chr; long w_start, w_end;
    int nL, nR, nC, nM;
    std::vector<GenoMat> X;        // all components (alpha-weighted, tr=n)
    std::vector<int> sig;          // indices in X that are signal (middle categories)
    std::vector<int> sig_cat;      // category id (0..n_cat, n_cat = uncategorized)
    std::vector<int> sig_m;        // #SNPs per signal component
    int mid_total;                 // total middle SNPs actually partitioned
};

// One window's RAW cells (standardized, but NOT yet alpha-weighted / trace-
// normalized / partitioned). Produced serially by the main thread (only
// readBedBlock + standardization are serial); the heavy per-window building is
// done later in the parallel worker via build_part_window().
struct RawWindow {
    std::string chr; long w_start, w_end;
    int nL, nM, nR, nC;
    Cell wesL, wesM, wesR;         // standardized WES cells (left/middle/right)
    Cell com;                      // standardized, concatenated common cell (may be empty)
    int mid_i0;                    // .bim index of the middle cell's first SNP
};

// Rolling reader: partitions each chromosome into variable-width CELLS (>= W bp
// AND >= min_snps SNPs; see build_cells), caches L/M/R WES cells and L/M/R
// common cells, shifts them, and hands each window its cells as a RawWindow --
// WITHOUT building components (that happens in the worker). Because middle and
// flanks are all the SAME cells, each cell is read exactly once as it rolls
// right -> middle -> left.
class PartStream {
public:
    PartStream(const RunContext& ctx, long W, int min_snps)
        : ctx_(ctx), W_(W), min_snps_(min_snps), ci_(0), have_(false) {}

    bool next(RawWindow& rw) {
        while (true) {
            if (!have_) { if (ci_ >= ctx_.chr_order.size()) return false; setup_chrom(); have_ = true; }
            if (g_ > glast_) { ++ci_; have_ = false; continue; }
            bool produced = false;
            if (wesM_.X.cols() > 0) { assemble(rw); produced = true; }
            shift();
            if (produced) return true;
        }
    }
private:
    const RunContext& ctx_; long W_; int min_snps_; size_t ci_; bool have_;
    int c_lo_, c_hi_, cc_lo_, cc_hi_; long g_, glast_;
    std::vector< std::pair<int,int> > cells_;    // WES SNP-index ranges for this chr
    Cell wesL_, wesM_, wesR_, comL_, comM_, comR_;

    Cell empty_cell() const { Cell c; c.X = GenoMat(ctx_.n_inds, 0); c.i0 = -1; return c; }

    Cell wesCell(long p) const {
        if (p < 0 || p >= (long) cells_.size()) return empty_cell();
        return read_cell_idx(ctx_.wes_prefix, ctx_.wes_n_total, ctx_.wes_n_snps,
                             ctx_.geno_keep, cells_[p].first, cells_[p].second + 1);
    }
    Cell comCell(long p) const {
        if (!ctx_.use_common || p < 0 || p >= (long) cells_.size()) return empty_cell();
        long start = ctx_.wes_bim.bp[cells_[p].first];
        long next  = (p + 1 < (long) cells_.size()) ? ctx_.wes_bim.bp[cells_[p + 1].first]
                                                     : (ctx_.wes_bim.bp[c_hi_] + 1);
        return read_cell(ctx_.common_prefix, ctx_.common_n_total, ctx_.common_n_snps,
                         ctx_.common_bim, ctx_.common_keep, cc_lo_, cc_hi_, start, next - start);
    }
    void setup_chrom() {
        const std::string& chr = ctx_.chr_order[ci_];
        c_lo_ = ctx_.chr_lo[ci_]; c_hi_ = ctx_.chr_hi[ci_];
        common_chr_range(ctx_, chr, cc_lo_, cc_hi_);
        cells_ = build_cells(ctx_.wes_bim.bp, c_lo_, c_hi_, W_, min_snps_);
        g_ = 0; glast_ = (long) cells_.size() - 1;
        wesL_ = wesCell(-1); wesM_ = wesCell(0); wesR_ = wesCell(1);
        comL_ = comCell(-1); comM_ = comCell(0); comR_ = comCell(1);
    }
    void shift() {
        ++g_;
        if (g_ > glast_) return;
        wesL_ = std::move(wesM_); wesM_ = std::move(wesR_); wesR_ = wesCell(g_ + 1);
        comL_ = std::move(comM_); comM_ = std::move(comR_); comR_ = comCell(g_ + 1);
    }

    void assemble(RawWindow& rw) {
        rw = RawWindow();
        rw.chr = ctx_.chr_order[ci_];
        rw.w_start = ctx_.wes_bim.bp[cells_[g_].first];
        rw.w_end   = ctx_.wes_bim.bp[cells_[g_].second];
        rw.mid_i0 = wesM_.i0;
        rw.nL = (int) wesL_.X.cols();
        rw.nM = (int) wesM_.X.cols();
        rw.nR = (int) wesR_.X.cols();
        rw.wesL = wesL_; rw.wesM = wesM_; rw.wesR = wesR_;    // copies (standardized)

        // concatenate the 3 common cells into ONE standardized cell (single copy;
        // the alpha weighting / trace-norm happen later, in the worker)
        int cm = (int) comL_.X.cols() + (int) comM_.X.cols() + (int) comR_.X.cols();
        rw.nC = cm; rw.com.i0 = -1;
        if (ctx_.use_common && cm > 0) {
            rw.com.X = GenoMat(ctx_.n_inds, cm);
            rw.com.maf.resize(cm);
            int off = 0; const Cell* cc[3] = { &comL_, &comM_, &comR_ };
            for (int t = 0; t < 3; ++t) {
                int mt = (int) cc[t]->X.cols();
                if (mt > 0) { rw.com.X.middleCols(off, mt) = cc[t]->X;
                              for (int j = 0; j < mt; ++j) rw.com.maf[off + j] = cc[t]->maf[j];
                              off += mt; }
            }
        } else {
            rw.com.X = GenoMat(ctx_.n_inds, 0);
        }
    }
};

// Build the component matrices for one window from its raw cells. Runs in the
// WORKER thread (thread-safe: only reads ctx.snp_cats/alpha/n_inds and calls the
// pure build_component). This is where the alpha weighting, trace-normalization
// and category partitioning happen -- moved off the serial main thread so they
// parallelize across cores.
static void build_part_window(const RawWindow& rw, const RunContext& ctx, PartWindow& pw) {
    pw = PartWindow();
    pw.chr = rw.chr; pw.w_start = rw.w_start; pw.w_end = rw.w_end;
    pw.nL = rw.nL; pw.nR = rw.nR; pw.nC = rw.nC; pw.nM = rw.nM;

    // per-SNP LD weights for the WES components (NULL if none were supplied)
    const std::vector<float>* wp = ctx.snp_weights.empty() ? 0 : &ctx.snp_weights;

    // background components: left, right (WES flanks), common
    if (rw.nL > 0) { std::vector<int> cols(rw.nL); for (int j = 0; j < rw.nL; ++j) cols[j] = j;
        GenoMat g = build_component(rw.wesL.X, rw.wesL.maf, cols, ctx.alpha, ctx.n_inds, wp, rw.wesL.i0);
        if (g.cols() > 0) pw.X.push_back(std::move(g)); }
    if (rw.nR > 0) { std::vector<int> cols(rw.nR); for (int j = 0; j < rw.nR; ++j) cols[j] = j;
        GenoMat g = build_component(rw.wesR.X, rw.wesR.maf, cols, ctx.alpha, ctx.n_inds, wp, rw.wesR.i0);
        if (g.cols() > 0) pw.X.push_back(std::move(g)); }
    if (ctx.use_common && rw.nC > 0) { std::vector<int> cols(rw.nC); for (int j = 0; j < rw.nC; ++j) cols[j] = j;
        // common fileset has its own .bim, so the WES-aligned weights don't apply
        GenoMat g = build_component(rw.com.X, rw.com.maf, cols, ctx.alpha, ctx.n_inds);
        if (g.cols() > 0) pw.X.push_back(std::move(g)); }

    // middle: partition columns by category
    int mid_i0 = rw.mid_i0;
    int ncat_snp = (int) ctx.snp_cats.size();
    std::vector< std::vector<int> > cat_cols(ctx.n_cat + 1);   // +1 = uncategorized
    for (int j = 0; j < rw.nM; ++j) {
        int snp = mid_i0 + j;
        if (snp < 0 || snp >= ncat_snp) continue;              // guard (no stop() in worker)
        const std::vector<int>& cs = ctx.snp_cats[snp];
        if (cs.empty()) cat_cols[ctx.n_cat].push_back(j);
        else for (size_t k = 0; k < cs.size(); ++k) cat_cols[cs[k]].push_back(j);
    }
    pw.mid_total = 0;
    for (int c = 0; c <= ctx.n_cat; ++c) {
        int m = (int) cat_cols[c].size();
        if (m == 0) continue;
        GenoMat g = build_component(rw.wesM.X, rw.wesM.maf, cat_cols[c], ctx.alpha, ctx.n_inds, wp, rw.wesM.i0);
        if (g.cols() == 0) continue;
        pw.sig.push_back((int) pw.X.size());
        pw.sig_cat.push_back(c);
        pw.sig_m.push_back(m);
        pw.mid_total += m;
        pw.X.push_back(std::move(g));
    }
}

// ===========================================================================
// Low-rank AI-REML kernels (copied from he_sliding_window.cpp; work on any set
// of component matrices via the concatenated Gram).
// ===========================================================================
struct RemlWindow {
    Eigen::MatrixXd Z, Gram; Eigen::VectorXd Zt1;
    std::vector<int> off, m; int C, K;
};
struct RemlEval { double logL; Eigen::VectorXd score; Eigen::MatrixXd AI; bool ok; };

static Eigen::VectorXd reml_wsqrt(const RemlWindow& W, const Eigen::VectorXd& s) {
    Eigen::VectorXd Ws(W.K);
    for (int l = 0; l < W.C; ++l) { double w = s[l] / W.m[l]; double sq = std::sqrt(w > 0 ? w : 0.0);
        for (int k = 0; k < W.m[l]; ++k) Ws[W.off[l] + k] = sq; }
    return Ws;
}
static Eigen::VectorXd reml_mom_start(const RemlWindow& W, const Eigen::VectorXd& Zty, double yty, int n) {
    int C = W.C, p = C + 1;
    Eigen::MatrixXd T = Eigen::MatrixXd::Zero(p, p); Eigen::VectorXd q(p);
    for (int a = 0; a < C; ++a) {
        for (int b = a; b < C; ++b) {
            double sfn = W.Gram.block(W.off[a], W.off[b], W.m[a], W.m[b]).squaredNorm();
            T(a, b) = sfn / ((double) W.m[a] * (double) W.m[b]); T(b, a) = T(a, b);
        }
        double tr = W.Gram.block(W.off[a], W.off[a], W.m[a], W.m[a]).trace() / W.m[a];
        T(a, C) = tr; T(C, a) = tr;
        q[a] = Zty.segment(W.off[a], W.m[a]).squaredNorm() / W.m[a];
    }
    T(C, C) = (double) n; q[C] = yty;
    return T.completeOrthogonalDecomposition().solve(q);
}
static double hutchinson_trace_Sinv(const Eigen::LLT<Eigen::MatrixXd>& llt, int K, int B, unsigned seed) {
    std::mt19937 rng(seed); std::uniform_int_distribution<int> bit(0, 1);
    double tr = 0.0; Eigen::VectorXd z(K);
    for (int b = 0; b < B; ++b) { for (int k = 0; k < K; ++k) z[k] = bit(rng) ? 1.0 : -1.0; tr += z.dot(llt.solve(z)); }
    return tr / B;
}
static bool reml_logL(const RemlWindow& W, const Eigen::VectorXd& Zty, const Eigen::VectorXd& y,
                      const Eigen::VectorXd& s, double& logL) {
    double se = s[W.C]; if (se <= 0) return false;
    Eigen::VectorXd Ws = reml_wsqrt(W, s);
    Eigen::MatrixXd S = (Ws.asDiagonal() * W.Gram * Ws.asDiagonal()) / se; S.diagonal().array() += 1.0;
    Eigen::LLT<Eigen::MatrixXd> llt(S); if (llt.info() != Eigen::Success) return false;
    double logdetS = 2.0 * llt.matrixLLT().diagonal().array().log().sum();
    Eigen::VectorXd g1 = Ws.cwiseProduct(W.Zt1), gy = Ws.cwiseProduct(Zty);
    Eigen::VectorXd v1 = Ws.cwiseProduct(llt.solve(g1)), vy = Ws.cwiseProduct(llt.solve(gy));
    Eigen::VectorXd Vi1 = Eigen::VectorXd::Ones(W.Z.rows()) / se - (W.Z * v1) / (se * se);
    Eigen::VectorXd Viy = y / se - (W.Z * vy) / (se * se);
    double c = Vi1.sum();
    Eigen::VectorXd r = Viy - Vi1 * (Viy.sum() / c);
    logL = -0.5 * ((double) W.Z.rows() * std::log(se) + logdetS + std::log(c) + y.dot(r));
    return true;
}
static RemlEval reml_eval(const RemlWindow& W, const Eigen::VectorXd& Zty, const Eigen::VectorXd& y,
                          const Eigen::VectorXd& s, unsigned seed) {
    RemlEval ev; ev.ok = false; int C = W.C, K = W.K, p = C + 1; double se = s[C];
    if (se <= 0) return ev; double Nn = (double) W.Z.rows();
    Eigen::VectorXd Ws = reml_wsqrt(W, s);
    Eigen::MatrixXd S = (Ws.asDiagonal() * W.Gram * Ws.asDiagonal()) / se; S.diagonal().array() += 1.0;
    Eigen::LLT<Eigen::MatrixXd> llt(S); if (llt.info() != Eigen::Success) return ev;
    double logdetS = 2.0 * llt.matrixLLT().diagonal().array().log().sum();
    Eigen::MatrixXd Gw = W.Gram * Ws.asDiagonal();
    Eigen::MatrixXd R = llt.solve(Gw.transpose());
    Eigen::MatrixXd ZVZ = W.Gram / se - (Gw * R) / (se * se);
    Eigen::VectorXd g1 = Ws.cwiseProduct(W.Zt1), gy = Ws.cwiseProduct(Zty);
    Eigen::VectorXd s1 = llt.solve(g1), sy = llt.solve(gy);
    Eigen::VectorXd ZV1 = W.Zt1 / se - (Gw * s1) / (se * se);
    Eigen::VectorXd Vi1 = Eigen::VectorXd::Ones((int) Nn) / se - (W.Z * Ws.cwiseProduct(s1)) / (se * se);
    Eigen::VectorXd Viy = y / se - (W.Z * Ws.cwiseProduct(sy)) / (se * se);
    double c = Vi1.sum();
    Eigen::VectorXd r = Viy - Vi1 * (Viy.sum() / c);
    Eigen::VectorXd Zr = W.Z.transpose() * r;
    Eigen::VectorXd a_r = Ws.cwiseProduct(Zr);
    Eigen::VectorXd Vir = r / se - (W.Z * Ws.cwiseProduct(llt.solve(a_r))) / (se * se);
    Eigen::VectorXd Pr = Vir - Vi1 * (Vir.sum() / c);
    Eigen::VectorXd ZPr = W.Z.transpose() * Pr;
    double trSinv = hutchinson_trace_Sinv(llt, K, 16, seed);
    double trVinv = Nn / se - (K - trSinv) / se;
    double trP = trVinv - Vi1.dot(Vi1) / c;
    Eigen::MatrixXd Q = ZVZ - (ZV1 * ZV1.transpose()) / c;
    Eigen::VectorXd score(p);
    for (int l = 0; l < C; ++l) {
        double trVG = 0, xv1 = 0, pl2 = 0;
        for (int k = 0; k < W.m[l]; ++k) { trVG += ZVZ(W.off[l] + k, W.off[l] + k);
            double z1 = ZV1[W.off[l] + k]; xv1 += z1 * z1; double zr = Zr[W.off[l] + k]; pl2 += zr * zr; }
        trVG /= W.m[l]; xv1 /= W.m[l]; pl2 /= W.m[l];
        score[l] = -0.5 * ((trVG - xv1 / c) - pl2);
    }
    score[C] = -0.5 * (trP - r.dot(r));
    Eigen::MatrixXd AI(p, p);
    for (int a = 0; a < C; ++a) for (int b = a; b < C; ++b) {
        Eigen::VectorXd pa = Zr.segment(W.off[a], W.m[a]), pb = Zr.segment(W.off[b], W.m[b]);
        Eigen::MatrixXd Qs = Q.block(W.off[a], W.off[b], W.m[a], W.m[b]);
        double v = 0.5 / ((double) W.m[a] * W.m[b]) * pa.dot(Qs * pb); AI(a, b) = v; AI(b, a) = v;
    }
    for (int l = 0; l < C; ++l) {
        Eigen::VectorXd pl = Zr.segment(W.off[l], W.m[l]), zp = ZPr.segment(W.off[l], W.m[l]);
        double v = 0.5 / W.m[l] * pl.dot(zp); AI(l, C) = v; AI(C, l) = v;
    }
    AI(C, C) = 0.5 * r.dot(Pr);
    ev.logL = -0.5 * (Nn * std::log(se) + logdetS + std::log(c) + y.dot(r));
    ev.score = score; ev.AI = AI; ev.ok = true; return ev;
}
static void reml_fit(const RemlWindow& W, const Eigen::VectorXd& Zty, const Eigen::VectorXd& y,
                     int max_iter, double tol, unsigned seed,
                     Eigen::VectorXd& s, bool& converged, int& iters, Eigen::MatrixXd& AI_out) {
    int p = W.C + 1;
    double vy = (y.array() - y.mean()).square().sum() / (y.size() - 1);
    double floor_v = 1e-6 * vy;
    double yty = y.squaredNorm();
    s = reml_mom_start(W, Zty, yty, (int) y.size());
    for (int l = 0; l < p; ++l) if (s[l] < floor_v) s[l] = floor_v;
    converged = false; iters = 0; AI_out = Eigen::MatrixXd::Zero(p, p);
    RemlEval cur = reml_eval(W, Zty, y, s, seed);
    if (!cur.ok) { s.setConstant(vy / p); s[W.C] = vy; cur = reml_eval(W, Zty, y, s, seed); if (!cur.ok) return; }
    AI_out = cur.AI;
    for (int it = 1; it <= max_iter; ++it) {
        iters = it;
        Eigen::VectorXd delta = cur.AI.ldlt().solve(cur.score);
        double step = 1.0; bool acc = false; Eigen::VectorXd s_acc;
        for (int h = 0; h < 25; ++h) {
            Eigen::VectorXd st = s + step * delta;
            for (int l = 0; l < p; ++l) if (st[l] < floor_v) st[l] = floor_v;
            double ll; if (reml_logL(W, Zty, y, st, ll) && ll >= cur.logL - 1e-8) { s_acc = st; acc = true; break; }
            step *= 0.5;
        }
        if (!acc) break;
        s = s_acc;
        RemlEval full = reml_eval(W, Zty, y, s, seed + (unsigned) it * 2654435761u);
        if (!full.ok) break;
        AI_out = full.AI;
        double dl = full.logL - cur.logL; cur = full;
        if (std::abs(dl) < tol) { converged = true; break; }
    }
}

// ===========================================================================
// Per-window results + compute kernels
// ===========================================================================
struct PartResult {
    // per signal component (category), per trait
    std::vector< std::vector<double> > vg, se_vg, h2, se_h2;   // [sig][trait]
    std::vector< std::vector<int> > conv, iters;               // REML
    // per trait (whole-window model fit): HE = genetic part of sigma'q
    // (sum_c sigma_c q_c, residual term excluded; higher = better),
    // REML = converged restricted log-likelihood.
    std::vector<double> fit;
    // window metadata (copied from the built PartWindow so the main-thread
    // driver can write output without re-touching the components)
    std::string chr; long w_start, w_end; int nL, nR, nC, mid_total;
    std::vector<int> sig_cat, sig_m;
};

static void alloc_result(PartResult& res, const PartWindow& pw, int P) {
    int nsig = (int) pw.sig.size();
    res.vg.assign(nsig, std::vector<double>(P, NA_REAL));
    res.se_vg.assign(nsig, std::vector<double>(P, NA_REAL));
    res.h2.assign(nsig, std::vector<double>(P, NA_REAL));
    res.se_h2.assign(nsig, std::vector<double>(P, NA_REAL));
    res.conv.assign(nsig, std::vector<int>(P, 0));
    res.iters.assign(nsig, std::vector<int>(P, 0));
    res.fit.assign(P, NA_REAL);
    res.chr = pw.chr; res.w_start = pw.w_start; res.w_end = pw.w_end;
    res.nL = pw.nL; res.nR = pw.nR; res.nC = pw.nC; res.mid_total = pw.mid_total;
    res.sig_cat = pw.sig_cat; res.sig_m = pw.sig_m;
}

// ---- HE (randomized MoM), records every signal component ----
static void he_part_compute(const PartWindow& pw, const Eigen::MatrixXd& Y,
                            const std::vector<double>& Vp, int nmcmc, unsigned seed,
                            bool se, PartResult& res) {
    int C = (int) pw.X.size(), env = C, n = (int) Y.rows(), P = (int) Y.cols();
    int nsig = (int) pw.sig.size();
    alloc_result(res, pw, P);
    if (C == 0) return;

    std::vector<double> M(C);
    for (int c = 0; c < C; ++c) M[c] = (double) pw.X[c].cols();

    GenoMat Yf = Y.cast<float>();
    std::vector<GenoMat> XtY(C);
    Eigen::MatrixXd T = Eigen::MatrixXd::Zero(C + 1, C + 1);
    for (int c = 0; c < C; ++c) {
        XtY[c] = pw.X[c].transpose() * Yf;
        double tr = pw.X[c].colwise().squaredNorm().cast<double>().sum() / M[c];   // tr(K_c) (~n, exact)
        T(c, env) = tr; T(env, c) = tr;
    }
    T(env, env) = n;

    int B = nmcmc;
    GenoMat Z(n, B);
    std::mt19937 rng(seed); std::normal_distribution<float> nd(0.f, 1.f);
    for (int p = 0; p < B; ++p) for (int i = 0; i < n; ++i) Z(i, p) = nd(rng);
    std::vector<GenoMat> Gz(C);
    for (int c = 0; c < C; ++c) Gz[c] = (pw.X[c] * (pw.X[c].transpose() * Z)) / (float) M[c];
    for (int a = 0; a < C; ++a) for (int b = a; b < C; ++b) {
        double t = 0.0; for (int p = 0; p < B; ++p) t += (double) Gz[a].col(p).dot(Gz[b].col(p));
        T(a, b) = t / B; T(b, a) = T(a, b);
    }
    Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> Tcod(T);
    Eigen::MatrixXd Tinv = Tcod.pseudoInverse();

    for (int t = 0; t < P; ++t) {
        Eigen::VectorXd q(C + 1);
        for (int c = 0; c < C; ++c) q[c] = XtY[c].col(t).cast<double>().squaredNorm() / M[c];
        q[env] = Y.col(t).squaredNorm();
        Eigen::VectorXd sigma = Tcod.solve(q);
        // Model fit for choosing alpha: the GENETIC part of sigma'q only.
        // The full sigma'q would add sigma[env] * y'y, and since y'y is identical
        // for every alpha and huge (~n per window, summed over thousands of
        // windows) it swamps the genetic signal -- maximizing it then just
        // rewards whichever alpha leaves the most variance in the residual, which
        // biases the selected alpha upward (validated in simulation). Dropping the
        // residual term restores calibration.
        double fit_gen = 0.0;
        for (int c = 0; c < C; ++c) fit_gen += sigma[c] * q[c];
        res.fit[t] = fit_gen;

        Eigen::MatrixXd Cov;
        if (se) {
            GenoMat u = (float) 0.0 * Z;                     // will hold Sigma z
            {   // apply Sigma = sum_c sigma[c]/M[c] X_c X_c' + sigma[env] I
                u = (float) sigma[env] * Z;
                for (int c = 0; c < C; ++c) u.noalias() += (float)(sigma[c] / M[c]) * (pw.X[c] * (pw.X[c].transpose() * Z));
            }
            std::vector<GenoMat> skju(C + 1);
            for (int l = 0; l <= C; ++l) {
                GenoMat Klu;
                if (l == C) Klu = u; else Klu = (pw.X[l] * (pw.X[l].transpose() * u)) / (float) M[l];
                GenoMat sk = (float) sigma[env] * Klu;
                for (int c = 0; c < C; ++c) sk.noalias() += (float)(sigma[c] / M[c]) * (pw.X[c] * (pw.X[c].transpose() * Klu));
                skju[l] = std::move(sk);
            }
            Eigen::MatrixXd Covq(C + 1, C + 1);
            for (int i = 0; i <= C; ++i) {
                const GenoMat& kizi = (i < C) ? Gz[i] : Z;
                for (int j = 0; j <= C; ++j) {
                    double v = 0.0; for (int p = 0; p < B; ++p) v += (double) kizi.col(p).dot(skju[j].col(p));
                    Covq(i, j) = v;
                }
            }
            Covq = (Covq / B) * 2.0; Covq = (0.5 * (Covq + Covq.transpose())).eval();
            Cov = Tinv * Covq * Tinv;
        }

        for (int k = 0; k < nsig; ++k) {
            int c = pw.sig[k];
            double vg = sigma[c];                            // = h2 contribution (tr(K)=n)
            res.vg[k][t] = vg;
            res.h2[k][t] = (Vp[t] > 0) ? vg / Vp[t] : NA_REAL;
            if (se) {
                if (Cov(c, c) > 0) {
                    double sev = std::sqrt(Cov(c, c));
                    res.se_vg[k][t] = sev;
                    res.se_h2[k][t] = (Vp[t] > 0) ? sev / Vp[t] : NA_REAL;
                }
            }
        }
    }
}

// ---- REML (low-rank), records every signal component ----
static void reml_part_compute(const PartWindow& pw, const Eigen::MatrixXd& Y,
                              const std::vector<double>& Vp, int max_iter, double tol,
                              bool se, unsigned seed, PartResult& res) {
    int C = (int) pw.X.size(), n = (int) Y.rows(), P = (int) Y.cols();
    int nsig = (int) pw.sig.size();
    alloc_result(res, pw, P);
    if (C == 0) return;

    // Build concatenated RemlWindow (Gram via float rankUpdate -> double)
    RemlWindow rw; rw.C = C; rw.off.resize(C); rw.m.resize(C);
    int K = 0; for (int c = 0; c < C; ++c) { rw.off[c] = K; rw.m[c] = (int) pw.X[c].cols(); K += rw.m[c]; }
    rw.K = K;
    GenoMat Zf(n, K);
    for (int c = 0; c < C; ++c) Zf.middleCols(rw.off[c], rw.m[c]) = pw.X[c];
    GenoMat Gf = GenoMat::Zero(K, K);
    Gf.selfadjointView<Eigen::Upper>().rankUpdate(Zf.transpose());
    rw.Gram = GenoMat(Gf.selfadjointView<Eigen::Upper>()).cast<double>();
    rw.Z = Zf.cast<double>();
    rw.Zt1 = rw.Z.colwise().sum().transpose();

    for (int t = 0; t < P; ++t) {
        Eigen::VectorXd y = Y.col(t);
        Eigen::VectorXd Zty = rw.Z.transpose() * y;
        Eigen::VectorXd s; bool conv; int iters; Eigen::MatrixXd AI;
        reml_fit(rw, Zty, y, max_iter, tol, seed + (unsigned) t * 7919u, s, conv, iters, AI);
        if (conv) { double ll; if (reml_logL(rw, Zty, y, s, ll)) res.fit[t] = ll; }  // restricted logL at solution
        Eigen::MatrixXd Cov;
        if (se && conv) Cov = AI.completeOrthogonalDecomposition().pseudoInverse();
        for (int k = 0; k < nsig; ++k) {
            int c = pw.sig[k];
            double vg = s[c];
            res.vg[k][t] = vg;
            res.h2[k][t] = (conv && Vp[t] > 0) ? vg / Vp[t] : NA_REAL;
            res.conv[k][t] = conv ? 1 : 0; res.iters[k][t] = iters;
            if (se && conv && Cov(c, c) > 0) {
                double sev = std::sqrt(Cov(c, c));
                res.se_vg[k][t] = sev;
                res.se_h2[k][t] = (Vp[t] > 0) ? sev / Vp[t] : NA_REAL;
            }
        }
    }
}

// ===========================================================================
// Progress bar (same throttled single-line style as he_sliding_window.cpp)
// ===========================================================================
// Count windows (= cells) that will be processed, without genotype I/O, so we
// can show a percentage / ETA. Uses the same variable-cell partition as the run.
static int count_total_windows(const RunContext& ctx, long W, int min_snps) {
    int total = 0;
    for (size_t ci = 0; ci < ctx.chr_order.size(); ++ci)
        total += (int) build_cells(ctx.wes_bim.bp, ctx.chr_lo[ci], ctx.chr_hi[ci], W, min_snps).size();
    return total;
}

// mm:ss formatter for elapsed / remaining time.
static std::string fmt_time(double s) {
    if (s < 0 || !(s == s)) s = 0;               // guard NaN / negative
    int t = (int)(s + 0.5);
    char b[24];
    std::snprintf(b, sizeof(b), "%d:%02d", t / 60, t % 60);
    return std::string(b);
}

// Throttled, single-line progress bar with elapsed time and ETA.
struct Progress {
    typedef std::chrono::steady_clock clock;
    int total, done;
    clock::time_point t0, tlast;
    std::string tag;

    void start(const std::string& t, int tot) {
        tag = t; total = tot; done = 0;
        t0 = tlast = clock::now();
        Rcout << "[" << tag << "] scanning " << total << " windows...\n";
    }
    static double secs(clock::time_point a, clock::time_point b) {
        return std::chrono::duration<double>(b - a).count();
    }
    void tick(const std::string& chr, int mid_snps) {
        ++done;
        clock::time_point now = clock::now();
        bool last = (done >= total);
        if (!last && secs(tlast, now) < 0.2) return;   // throttle to ~5 Hz
        tlast = now;

        double elapsed = secs(t0, now);
        double frac = total > 0 ? (double) done / total : 1.0;
        double eta = (frac > 0.0) ? elapsed * (1.0 - frac) / frac : 0.0;

        const int barw = 24;
        int filled = (int)(barw * frac + 0.5);
        std::string bar(barw, '-');
        for (int i = 0; i < filled && i < barw; ++i) bar[i] = '#';

        char buf[256];
        std::snprintf(buf, sizeof(buf),
            "\r[%s] [%s] %3d%% | chr %-3s | %d/%d | mid=%d SNPs | %s elapsed | ~%s left    ",
            tag.c_str(), bar.c_str(), (int)(100.0 * frac + 0.5),
            chr.c_str(), done, total, mid_snps,
            fmt_time(elapsed).c_str(), fmt_time(eta).c_str());
        Rcout << buf << std::flush;
        if (last) Rcout << "\n";
    }
};

// ===========================================================================
// Driver (parallel over windows in batches; streams a long-format table)
// ===========================================================================
struct PartParams {
    long W; int min_snps; int method;   // method: 0 = HE, 1 = REML
    int nmcmc; int max_iter; double tol; bool se;
    std::string out_file; int batch_size; int n_threads; unsigned seed;
};

// Workers now BUILD each window's components (build_part_window) and THEN
// estimate -- so the alpha weighting / trace-normalization / partitioning run
// in parallel across cores instead of serially on the main thread.
struct HEWorker : public RcppParallel::Worker {
    const std::vector<RawWindow>& batch; const RunContext& ctx;
    const Eigen::MatrixXd& Y; const std::vector<double>& Vp;
    const PartParams& pr; unsigned soff; std::vector<PartResult>& out;
    HEWorker(const std::vector<RawWindow>& b, const RunContext& ctx, const Eigen::MatrixXd& Y,
             const std::vector<double>& Vp, const PartParams& pr, unsigned soff, std::vector<PartResult>& out)
        : batch(b), ctx(ctx), Y(Y), Vp(Vp), pr(pr), soff(soff), out(out) {}
    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t w = begin; w < end; ++w) {
            PartWindow pw; build_part_window(batch[w], ctx, pw);
            he_part_compute(pw, Y, Vp, pr.nmcmc, pr.seed + soff + (unsigned) w, pr.se, out[w]);
        }
    }
};
struct RMLWorker : public RcppParallel::Worker {
    const std::vector<RawWindow>& batch; const RunContext& ctx;
    const Eigen::MatrixXd& Y; const std::vector<double>& Vp;
    const PartParams& pr; unsigned soff; std::vector<PartResult>& out;
    RMLWorker(const std::vector<RawWindow>& b, const RunContext& ctx, const Eigen::MatrixXd& Y,
              const std::vector<double>& Vp, const PartParams& pr, unsigned soff, std::vector<PartResult>& out)
        : batch(b), ctx(ctx), Y(Y), Vp(Vp), pr(pr), soff(soff), out(out) {}
    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t w = begin; w < end; ++w) {
            PartWindow pw; build_part_window(batch[w], ctx, pw);
            reml_part_compute(pw, Y, Vp, pr.max_iter, pr.tol, pr.se, pr.seed + soff + (unsigned) w, out[w]);
        }
    }
};

static Rcpp::List part_driver(RunContext& ctx, const PartParams& pr) {
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
    for (int t = 0; t < P; ++t) { double m = Y.col(t).mean(); Vp[t] = (Y.col(t).array() - m).square().sum() / (Y.rows() - 1); }

    std::ofstream fout;
    bool tofile = !pr.out_file.empty();
    auto wr = [](std::ofstream& f, double v) { if (std::isnan(v)) f << "NA"; else f << v; };
    if (tofile) {
        fout.open(pr.out_file.c_str());
        if (!fout.is_open()) stop("Could not open out_file: " + pr.out_file);
        fout << "chr\tstart\tend\tn_left\tn_right\tn_common\tphenotype\tcategory\tm_c\tvg\tse_vg\th2\tse_h2\tenrichment\tfit";
        if (pr.method == 1) fout << "\tconverged\tn_iter";
        fout << "\n";
    }
    // `fit` = whole-window model fit, repeated on each category row of a window:
    //   HE   -> GENETIC part of sigma'q, i.e. sum_c sigma_c q_c over the genetic
    //           components only (the residual term sigma_env * y'y is excluded --
    //           it is alpha-independent and would swamp the signal; see the note
    //           in he_part_compute). Larger = better.
    //   REML -> converged restricted log-likelihood (larger = better)
    // Summed over windows it is the per-trait score for choosing alpha.

    // genome-wide accumulators: per category (incl uncategorized) per trait
    int NC = ctx.n_cat + 1;
    std::vector< std::vector<double> > tot_vg(NC, std::vector<double>(P, 0.0));
    std::vector< std::vector<double> > tot_var(NC, std::vector<double>(P, 0.0));
    std::vector<double> tot_fit(P, 0.0);        // summed fit / logL per trait

    PartStream stream(ctx, pr.W, pr.min_snps);
    long processed = 0; int n_win = 0;
    Progress prog; prog.start(pr.method == 0 ? "HE-part" : "REML-part",
                              count_total_windows(ctx, pr.W, pr.min_snps));

    // -------- DEBUG timing/size checkpoints (set DEBUG_PART 1 to enable) ------
#define DEBUG_PART 0
    typedef std::chrono::steady_clock dbgclk;
    double t_read = 0.0, t_comp = 0.0, t_write = 0.0;   // seconds in each phase
    long long snps_seen = 0; int max_K = 0, max_nC = 0;  // component-size peaks
    int dbg_batch = 0;
    (void) t_read; (void) t_comp; (void) t_write;        // silence when DEBUG_PART=0
    (void) snps_seen; (void) max_K; (void) max_nC; (void) dbg_batch;
    auto dbg_secs = [](dbgclk::time_point a, dbgclk::time_point b) {
        return std::chrono::duration<double>(b - a).count(); };

    while (true) {
        dbgclk::time_point t0 = dbgclk::now();
        std::vector<RawWindow> batch;
        RawWindow rw;
        while ((int) batch.size() < pr.batch_size && stream.next(rw)) batch.push_back(std::move(rw));
        if (batch.empty()) break;
        dbgclk::time_point t1 = dbgclk::now();
        t_read += dbg_secs(t0, t1);   // read + standardize only (build now in worker)
#if DEBUG_PART
        for (size_t bb = 0; bb < batch.size(); ++bb) {
            int Kt = batch[bb].nL + batch[bb].nM + batch[bb].nR + batch[bb].nC;
            snps_seen += Kt; if (Kt > max_K) max_K = Kt; if (batch[bb].nC > max_nC) max_nC = batch[bb].nC;
        }
        if (dbg_batch < 3) {
            const RawWindow& w0 = batch.front();
            Rcout << "[dbg] batch " << dbg_batch << " (" << batch.size() << " win) first: "
                  << "chr " << w0.chr << " nL=" << w0.nL << " nM=" << w0.nM
                  << " nR=" << w0.nR << " nC=" << w0.nC
                  << " | read+std so far=" << t_read << "s\n";
        }
#endif

        std::vector<PartResult> results(batch.size());
        dbgclk::time_point t2 = dbgclk::now();
        if (pr.method == 0) { HEWorker wk(batch, ctx, Y, Vp, pr, (unsigned) processed, results);
                              RcppParallel::parallelFor(0, batch.size(), wk); }
        else                { RMLWorker wk(batch, ctx, Y, Vp, pr, (unsigned) processed, results);
                              RcppParallel::parallelFor(0, batch.size(), wk); }
        dbgclk::time_point t3 = dbgclk::now();
        t_comp += dbg_secs(t2, t3);   // build + estimate (parallel)
#if DEBUG_PART
        if (dbg_batch < 3 || (dbg_batch % 50 == 0))
            Rcout << "[dbg] cum: read+std=" << t_read << "s  build+compute=" << t_comp
                  << "s  write=" << t_write << "s  windows=" << n_win
                  << "  max_K=" << max_K << "  max_nC=" << max_nC << "\n";
        ++dbg_batch;
        dbgclk::time_point tw0 = dbgclk::now();
#endif

        for (size_t b = 0; b < batch.size(); ++b) {
            const PartResult& r = results[b];
            int nsig = (int) r.sig_cat.size();
            for (int t = 0; t < P; ++t) {
                std::string ph = as<std::string>(ctx.trait_names[t]);
                if (!std::isnan(r.fit[t])) tot_fit[t] += r.fit[t];   // once per window x trait
                // per-trait middle total h2 (sum of signal vg) for enrichment
                double mid_vg = 0.0;
                for (int k = 0; k < nsig; ++k) if (!std::isnan(r.vg[k][t])) mid_vg += r.vg[k][t];
                for (int k = 0; k < nsig; ++k) {
                    int cat = r.sig_cat[k]; int mc = r.sig_m[k];
                    double vg = r.vg[k][t];
                    double expected = (r.mid_total > 0) ? (double) mc / r.mid_total : NA_REAL;
                    double share = (mid_vg != 0.0 && !std::isnan(vg)) ? vg / mid_vg : NA_REAL;
                    double enr = (!std::isnan(share) && !std::isnan(expected) && expected > 0) ? share / expected : NA_REAL;
                    if (tofile) {
                        fout << r.chr << '\t' << r.w_start << '\t' << r.w_end << '\t'
                             << r.nL << '\t' << r.nR << '\t' << r.nC << '\t' << ph << '\t'
                             << ctx.cat_names[cat] << '\t' << mc << '\t';
                        wr(fout, vg); fout << '\t'; wr(fout, r.se_vg[k][t]); fout << '\t';
                        wr(fout, r.h2[k][t]); fout << '\t'; wr(fout, r.se_h2[k][t]); fout << '\t';
                        wr(fout, enr); fout << '\t'; wr(fout, r.fit[t]);
                        if (pr.method == 1) fout << '\t' << r.conv[k][t] << '\t' << r.iters[k][t];
                        fout << '\n';
                    }
                    if (!std::isnan(vg)) tot_vg[cat][t] += vg;
                    if (!std::isnan(r.se_vg[k][t])) tot_var[cat][t] += r.se_vg[k][t] * r.se_vg[k][t];
                }
            }
            ++n_win;
            prog.tick(r.chr, r.mid_total);          // one tick per finished window
        }
        processed += (long) batch.size();
        if (tofile) fout.flush();
        Rcpp::checkUserInterrupt();
#if DEBUG_PART
        t_write += dbg_secs(tw0, dbgclk::now());
#endif
    }

#if DEBUG_PART
    Rcout << "[dbg] TIMING TOTAL: read(+build)=" << t_read << "s  compute=" << t_comp
          << "s  write=" << t_write << "s\n"
          << "[dbg] windows=" << n_win << "  mean SNPs/window(K_total)="
          << (n_win > 0 ? (double) snps_seen / n_win : 0.0)
          << "  max_K=" << max_K << "  max_nC(common)=" << max_nC << "\n";
#endif

    // genome-wide TOTAL rows (one per category per trait)
    if (tofile) {
        for (int t = 0; t < P; ++t) {
            std::string ph = as<std::string>(ctx.trait_names[t]);
            for (int c = 0; c < NC; ++c) {
                double vg = tot_vg[c][t], sev = std::sqrt(tot_var[c][t]);
                fout << "TOTAL\tNA\tNA\tNA\tNA\tNA\t" << ph << '\t' << ctx.cat_names[c] << "\tNA\t";
                wr(fout, vg); fout << '\t'; wr(fout, sev); fout << '\t';
                wr(fout, Vp[t] > 0 ? vg / Vp[t] : NA_REAL); fout << '\t';
                wr(fout, Vp[t] > 0 ? sev / Vp[t] : NA_REAL); fout << "\tNA";     // enrichment
                fout << '\t'; wr(fout, tot_fit[t]);                              // genome-wide fit / logL
                if (pr.method == 1) fout << "\tNA\tNA";
                fout << '\n';
            }
        }
        fout.close();
    }
    Rcout << "Windows estimated: " << n_win << "\n";

    // return genome-wide category h2 matrix (category x trait) for convenience
    Rcpp::NumericMatrix gh2(NC, P), gse(NC, P);
    CharacterVector rn(NC);
    for (int c = 0; c < NC; ++c) rn[c] = ctx.cat_names[c];
    for (int t = 0; t < P; ++t)
        for (int c = 0; c < NC; ++c) {
            gh2(c, t) = (Vp[t] > 0) ? tot_vg[c][t] / Vp[t] : NA_REAL;
            gse(c, t) = (Vp[t] > 0) ? std::sqrt(tot_var[c][t]) / Vp[t] : NA_REAL;
        }
    rownames(gh2) = rn; colnames(gh2) = ctx.trait_names;
    rownames(gse) = rn; colnames(gse) = ctx.trait_names;

    // genome-wide model fit per trait (HE: summed sigma'q; REML: summed logL).
    // Run once per alpha and pick the alpha maximizing this.
    Rcpp::NumericVector fit_out(P);
    for (int t = 0; t < P; ++t) fit_out[t] = tot_fit[t];
    fit_out.names() = ctx.trait_names;

    return List::create(_["genome_h2"] = gh2, _["genome_se"] = gse,
                        _["fit"] = fit_out,
                        _["categories"] = wrap(ctx.cat_names), _["trait_names"] = ctx.trait_names,
                        _["alpha"] = ctx.alpha, _["method"] = (pr.method == 0 ? "HE" : "REML"),
                        _["n_windows"] = n_win);
}

}  // end anonymous namespace (internal-linkage helpers)

// ===========================================================================
// Exports (external linkage; call into the anonymous-namespace implementation)
// ===========================================================================
// [[Rcpp::export]]
Rcpp::List he_sliding_window_part(const std::string& filename,
                                  const SEXP pheno_mat,
                                  const IntegerMatrix snp_cat,
                                  const CharacterVector cat_names,
                                  double window_size = 1e6,
                                  int min_snps = 1000,
                                  double alpha = -1.0,
                                  Rcpp::Nullable<Rcpp::String> common_filename = R_NilValue,
                                  int nmcmc = 20,
                                  bool se = true,
                                  std::string out_file = "",
                                  int batch_size = 8,
                                  int n_threads = 0,
                                  int seed = 12345,
                                  Rcpp::Nullable<Rcpp::NumericVector> weights = R_NilValue) {
    if (window_size <= 0) stop("window_size must be positive");
    if (min_snps < 1) min_snps = 1;
    RunContext ctx = setup_context(filename, pheno_mat, snp_cat, cat_names, alpha, common_filename, weights);
    PartParams pr; pr.W = (long) window_size; pr.min_snps = min_snps; pr.method = 0;
    pr.nmcmc = nmcmc; pr.max_iter = 0; pr.tol = 0.0; pr.se = se;
    pr.out_file = out_file; pr.batch_size = batch_size; pr.n_threads = n_threads; pr.seed = (unsigned) seed;
    return part_driver(ctx, pr);
}

// [[Rcpp::export]]
Rcpp::List reml_sliding_window_part(const std::string& filename,
                                    const SEXP pheno_mat,
                                    const IntegerMatrix snp_cat,
                                    const CharacterVector cat_names,
                                    double window_size = 1e6,
                                    int min_snps = 1000,
                                    double alpha = -1.0,
                                    Rcpp::Nullable<Rcpp::String> common_filename = R_NilValue,
                                    int max_iter = 100,
                                    double tol = 1e-4,
                                    bool se = true,
                                    std::string out_file = "",
                                    int batch_size = 4,
                                    int n_threads = 0,
                                    int seed = 12345,
                                    Rcpp::Nullable<Rcpp::NumericVector> weights = R_NilValue) {
    if (window_size <= 0) stop("window_size must be positive");
    if (min_snps < 1) min_snps = 1;
    RunContext ctx = setup_context(filename, pheno_mat, snp_cat, cat_names, alpha, common_filename, weights);
    PartParams pr; pr.W = (long) window_size; pr.min_snps = min_snps; pr.method = 1;
    pr.nmcmc = 0; pr.max_iter = max_iter; pr.tol = tol; pr.se = se;
    pr.out_file = out_file; pr.batch_size = batch_size; pr.n_threads = n_threads; pr.seed = (unsigned) seed;
    return part_driver(ctx, pr);
}
