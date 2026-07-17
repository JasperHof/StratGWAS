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

using namespace Rcpp;

// ===========================================================================
// FLEX-style PARTITIONED sliding-window heritability (HE) + a GENOME-WIDE
// background GRM.
//
// This is he_sliding_window_part()'s HE estimator, unchanged, PLUS one extra
// background variance component: a genome-wide GRM built from the SAME common
// fileset (common_filename) that the local common GRM already uses. Its job is
// to soak up the DISTAL genetic variance (causals on other chromosomes) that
// local flanks cannot capture, which is what leaks into and accumulates across
// null windows.
//
// The genome-wide GRM is NEVER materialized (no n x n matrix). Instead, in ONE
// streaming pass over the common SNPs we precompute only the quantities the HE
// moment system needs, using a SHARED probe block Z (the same probes every
// window uses):
//     U    = K_gw Z            (n x B, a few MB)
//     q_gw = y' K_gw y         (one scalar per trait)
//     trK  = tr(K_gw) = n      (trace-normalized, exact)
//     trK2 = tr(K_gw^2)        (~ ||U||_F^2 / B)
// Then every window reuses these with zero extra I/O:
//     tr(K_gw K_c) ~ (1/B) sum_b (K_c z_b) . U_b        (K_c z_b already formed)
// so K_gw slots in as one more row/column of T and one more entry of q.
// After the pass, only U (n x B) and a handful of scalars are retained; the
// streamed genotype blocks are freed as we go.
//
// Same input arguments as he_sliding_window_part(). If common_filename is NULL
// the genome-wide component is simply absent and this reduces to the plain HE
// partitioned estimator.
//
// NOTE: shared global probes. Unlike he_sliding_window_part (which seeds probes
// per window), here ALL windows use one probe block Z so the precomputed U is
// reusable. Still fully reproducible from `seed`.
//
// SE note: the middle-category SE is computed exactly as before over the
// window components + residual, treating the genome-wide background as a KNOWN
// offset (its own sampling uncertainty, O(1/sqrt(M_gw)) with M_gw ~ all common
// SNPs, is negligible). The point estimate fully conditions on it.
// ===========================================================================

namespace {

typedef Eigen::MatrixXf GenoMat;

// --------------------------------------------------------------------------
// Shared low-level helpers (identical to he_sliding_window_part)
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

struct Cell {
    GenoMat X;
    std::vector<float> maf;
    int i0;
};

static void standardize_capture_maf(GenoMat& X, std::vector<float>& maf) {
    const int n = (int) X.rows();
    maf.assign(X.cols(), 0.f);
    for (int col = 0; col < X.cols(); ++col) {
        double sum = 0.0; int n_valid = 0;
        for (int i = 0; i < n; ++i) { float g = X(i, col); if (g != -1.f) { sum += g; ++n_valid; } }
        if (n_valid == 0) { X.col(col).setZero(); maf[col] = 0.f; continue; }
        double mean = sum / n_valid;
        double f = mean / 2.0; if (f > 0.5) f = 1.0 - f;
        maf[col] = (float) f;
        double sq = 0.0;
        for (int i = 0; i < n; ++i) { double g = (X(i, col) == -1.f) ? mean : (double) X(i, col); sq += (g - mean) * (g - mean); }
        double sd = std::sqrt(sq / (n - 1));
        if (sd > 1e-10) { for (int i = 0; i < n; ++i) { double g = (X(i, col) == -1.f) ? mean : (double) X(i, col); X(i, col) = (float)((g - mean) / sd); } }
        else X.col(col).setZero();
    }
}

static Cell read_cell(const std::string& prefix, int n_total, int n_snps,
                      const BimInfo& bim, const std::vector<int>& keep,
                      int lo, int hi, long cell_start, long W) {
    Cell cell; cell.i0 = -1;
    if (lo < 0 || hi < lo) { cell.X = GenoMat(keep.size(), 0); return cell; }
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

// Read a raw SNP-INDEX range [lo,hi] (across the whole fileset, not by bp),
// subset to the analysis individuals, standardize, and apply the alpha weight.
// Used only by the one-pass genome-wide background precompute. Low memory: the
// caller streams in modest index blocks, so only one block is alive at a time.
static GenoMat read_block_weighted(const std::string& prefix, int n_total, int n_snps,
                                   const std::vector<int>& keep, int lo, int hi, double alpha) {
    int m = hi - lo + 1, n_keep = (int) keep.size();
    GenoMat A(n_keep, m);
    IntegerMatrix block = readBedBlock(prefix + ".bed", n_total, n_snps, 0, n_total - 1, lo, hi);
    for (int i = 0; i < n_keep; ++i)
        for (int j = 0; j < m; ++j)
            A(i, j) = (float) block(keep[i], j);
    std::vector<float> maf;
    standardize_capture_maf(A, maf);
    if (alpha != -1.0) {
        double e = (1.0 + alpha) / 2.0;
        for (int j = 0; j < m; ++j) {
            double f = maf[j];
            double w = (f > 0.0 && f < 1.0) ? std::pow(2.0 * f * (1.0 - f), e) : 1.0;
            A.col(j) *= (float) w;
        }
    }
    return A;
}

static GenoMat build_component(const GenoMat& X, const std::vector<float>& maf,
                               const std::vector<int>& cols, double alpha, int n) {
    int m = (int) cols.size();
    if (m == 0) return GenoMat(X.rows(), 0);
    GenoMat W(X.rows(), m);
    bool unit = (alpha == -1.0);
    double e = (1.0 + alpha) / 2.0;
    for (int k = 0; k < m; ++k) {
        int j = cols[k];
        if (unit) { W.col(k) = X.col(j); continue; }
        double f = maf[j];
        double w = (f > 0.0 && f < 1.0) ? std::pow(2.0 * f * (1.0 - f), e) : 1.0;
        W.col(k) = X.col(j) * (float) w;
    }
    double ss = W.colwise().squaredNorm().cast<double>().sum();
    if (ss < 1e-12) return GenoMat(X.rows(), 0);
    double c = std::sqrt((double) n * m / ss);
    W *= (float) c;
    return W;
}

// --------------------------------------------------------------------------
// Context (identical to he_sliding_window_part)
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
    int n_cat;
    std::vector<std::string> cat_names;
    std::vector< std::vector<int> > snp_cats;
    double alpha;
};

static RunContext setup_context(const std::string& filename, const SEXP pheno_mat,
                                const IntegerMatrix& snp_cat, const CharacterVector& cat_names,
                                double alpha, Rcpp::Nullable<Rcpp::String> common_filename) {
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
    ctx.cat_names.push_back("uncategorized");
    ctx.snp_cats.resize(ctx.wes_n_snps);
    for (int j = 0; j < ctx.wes_n_snps; ++j)
        for (int c = 0; c < ctx.n_cat; ++c)
            if (snp_cat(j, c) == 1) ctx.snp_cats[j].push_back(c);

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

    // Mean-center each phenotype (q_env = y'y assumes centered y; genome-wide
    // and local GRM columns are already centered so their q are unaffected).
    for (int j = 0; j < ctx.n_pheno; ++j) {
        double m = ctx.Y.col(j).mean();
        ctx.Y.col(j).array() -= m;
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

static void common_chr_range(const RunContext& ctx, const std::string& chr, int& cc_lo, int& cc_hi) {
    cc_lo = cc_hi = -1;
    if (!ctx.use_common) return;
    for (int j = 0; j < ctx.common_n_snps; ++j)
        if (ctx.common_bim.chr[j] == chr) { if (cc_lo == -1) cc_lo = j; cc_hi = j; }
}

// --------------------------------------------------------------------------
// Genome-wide background GRM: one-pass precompute (shared probes), low memory.
// Only U (n x B) + scalars are kept; every genotype block is freed as read.
// --------------------------------------------------------------------------
struct GwBackground {
    bool active = false;
    GenoMat U;                 // n x B  = K_gw Z
    double trK = 0.0;          // tr(K_gw)  (= n after trace-normalization)
    double trK2 = 0.0;         // tr(K_gw^2) ~ ||U||_F^2 / B
    std::vector<double> q;     // per trait: y' K_gw y
    long M = 0;                // # common SNPs used
};

static GwBackground precompute_gw(const RunContext& ctx, const GenoMat& Z, int block) {
    GwBackground gw;
    if (!ctx.use_common) return gw;
    const int n = ctx.n_inds, B = (int) Z.cols(), P = ctx.n_pheno;
    const int total = ctx.common_n_snps;

    Eigen::MatrixXd Vacc = Eigen::MatrixXd::Zero(n, B);   // sum_k A_k (A_k' Z), double
    double sumA2 = 0.0;                                   // ||A||_F^2
    std::vector<double> sy2(P, 0.0);                      // per trait ||A'y||^2
    GenoMat Yf = ctx.Y.cast<float>();                     // n x P
    long M = 0;

    Rcout << "[HE-gw] precomputing genome-wide background from " << total
          << " common SNPs (one pass, " << B << " shared probes)...\n";
    int done = 0;
    for (int lo = 0; lo < total; lo += block) {
        int hi = std::min(lo + block, total) - 1;
        GenoMat A = read_block_weighted(ctx.common_prefix, ctx.common_n_total,
                                        ctx.common_n_snps, ctx.common_keep, lo, hi, ctx.alpha);
        sumA2 += A.squaredNorm();                          // double accumulate
        GenoMat AtZ = A.transpose() * Z;                   // m x B
        Vacc.noalias() += (A * AtZ).cast<double>();        // n x B
        for (int t = 0; t < P; ++t) {
            GenoMat Aty = A.transpose() * Yf.col(t);       // m x 1
            sy2[t] += (double) Aty.squaredNorm();
        }
        M += (long) A.cols();
        done = hi + 1;
        Rcpp::checkUserInterrupt();
        Rcout << "\r[HE-gw]   " << done << "/" << total << " common SNPs   " << std::flush;
    }
    Rcout << "\n";

    if (sumA2 < 1e-12) { gw.active = false; return gw; }
    double kappa = (double) n / sumA2;                     // K_gw = kappa * A A'  -> tr = n
    gw.U   = (kappa * Vacc).cast<float>();                 // K_gw Z
    gw.trK = (double) n;
    gw.trK2 = gw.U.cast<double>().squaredNorm() / B;        // ~ tr(K_gw^2)
    gw.q.assign(P, 0.0);
    for (int t = 0; t < P; ++t) gw.q[t] = kappa * sy2[t];  // y' K_gw y
    gw.M = M;
    gw.active = true;
    Rcout << "[HE-gw] genome-wide background ready (M=" << M
          << ", tr(K_gw)=" << gw.trK << ", tr(K_gw^2)~" << gw.trK2 << ")\n";
    return gw;
}

// --------------------------------------------------------------------------
// Window structs / stream (identical to he_sliding_window_part)
// --------------------------------------------------------------------------
struct PartWindow {
    std::string chr; long w_start;
    int nL, nR, nC, nM;
    std::vector<GenoMat> X;
    std::vector<int> sig;
    std::vector<int> sig_cat;
    std::vector<int> sig_m;
    int mid_total;
};

struct RawWindow {
    std::string chr; long w_start;
    int nL, nM, nR, nC;
    Cell wesL, wesM, wesR;
    Cell com;
    int mid_i0;
};

class PartStream {
public:
    PartStream(const RunContext& ctx, long W) : ctx_(ctx), W_(W), ci_(0), have_(false) {}
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
    const RunContext& ctx_; long W_; size_t ci_; bool have_;
    int c_lo_, c_hi_, cc_lo_, cc_hi_; long g_, glast_;
    Cell wesL_, wesM_, wesR_, comL_, comM_, comR_;

    Cell wesCell(long cell) const {
        return read_cell(ctx_.wes_prefix, ctx_.wes_n_total, ctx_.wes_n_snps,
                         ctx_.wes_bim, ctx_.geno_keep, c_lo_, c_hi_, cell * W_, W_);
    }
    Cell comCell(long cell) const {
        if (!ctx_.use_common) { Cell c; c.X = GenoMat(ctx_.n_inds, 0); c.i0 = -1; return c; }
        return read_cell(ctx_.common_prefix, ctx_.common_n_total, ctx_.common_n_snps,
                         ctx_.common_bim, ctx_.common_keep, cc_lo_, cc_hi_, cell * W_, W_);
    }
    void setup_chrom() {
        const std::string& chr = ctx_.chr_order[ci_];
        c_lo_ = ctx_.chr_lo[ci_]; c_hi_ = ctx_.chr_hi[ci_];
        common_chr_range(ctx_, chr, cc_lo_, cc_hi_);
        long min_bp = ctx_.wes_bim.bp[c_lo_], max_bp = ctx_.wes_bim.bp[c_hi_];
        g_ = min_bp / W_; glast_ = max_bp / W_;
        wesL_ = wesCell(g_ - 1); wesM_ = wesCell(g_); wesR_ = wesCell(g_ + 1);
        comL_ = comCell(g_ - 1); comM_ = comCell(g_); comR_ = comCell(g_ + 1);
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
        rw.w_start = g_ * W_;
        rw.mid_i0 = wesM_.i0;
        rw.nL = (int) wesL_.X.cols();
        rw.nM = (int) wesM_.X.cols();
        rw.nR = (int) wesR_.X.cols();
        rw.wesL = wesL_; rw.wesM = wesM_; rw.wesR = wesR_;
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

static void build_part_window(const RawWindow& rw, const RunContext& ctx, PartWindow& pw) {
    pw = PartWindow();
    pw.chr = rw.chr; pw.w_start = rw.w_start;
    pw.nL = rw.nL; pw.nR = rw.nR; pw.nC = rw.nC; pw.nM = rw.nM;

    if (rw.nL > 0) { std::vector<int> cols(rw.nL); for (int j = 0; j < rw.nL; ++j) cols[j] = j;
        GenoMat g = build_component(rw.wesL.X, rw.wesL.maf, cols, ctx.alpha, ctx.n_inds);
        if (g.cols() > 0) pw.X.push_back(std::move(g)); }
    if (rw.nR > 0) { std::vector<int> cols(rw.nR); for (int j = 0; j < rw.nR; ++j) cols[j] = j;
        GenoMat g = build_component(rw.wesR.X, rw.wesR.maf, cols, ctx.alpha, ctx.n_inds);
        if (g.cols() > 0) pw.X.push_back(std::move(g)); }
    if (ctx.use_common && rw.nC > 0) { std::vector<int> cols(rw.nC); for (int j = 0; j < rw.nC; ++j) cols[j] = j;
        GenoMat g = build_component(rw.com.X, rw.com.maf, cols, ctx.alpha, ctx.n_inds);
        if (g.cols() > 0) pw.X.push_back(std::move(g)); }

    int mid_i0 = rw.mid_i0;
    int ncat_snp = (int) ctx.snp_cats.size();
    std::vector< std::vector<int> > cat_cols(ctx.n_cat + 1);
    for (int j = 0; j < rw.nM; ++j) {
        int snp = mid_i0 + j;
        if (snp < 0 || snp >= ncat_snp) continue;
        const std::vector<int>& cs = ctx.snp_cats[snp];
        if (cs.empty()) cat_cols[ctx.n_cat].push_back(j);
        else for (size_t k = 0; k < cs.size(); ++k) cat_cols[cs[k]].push_back(j);
    }
    pw.mid_total = 0;
    for (int c = 0; c <= ctx.n_cat; ++c) {
        int m = (int) cat_cols[c].size();
        if (m == 0) continue;
        GenoMat g = build_component(rw.wesM.X, rw.wesM.maf, cat_cols[c], ctx.alpha, ctx.n_inds);
        if (g.cols() == 0) continue;
        pw.sig.push_back((int) pw.X.size());
        pw.sig_cat.push_back(c);
        pw.sig_m.push_back(m);
        pw.mid_total += m;
        pw.X.push_back(std::move(g));
    }
}

// --------------------------------------------------------------------------
// Per-window results
// --------------------------------------------------------------------------
struct PartResult {
    std::vector< std::vector<double> > vg, se_vg, h2, se_h2;
    std::vector<double> fit;
    std::string chr; long w_start; int nL, nR, nC, mid_total;
    std::vector<int> sig_cat, sig_m;
};

static void alloc_result(PartResult& res, const PartWindow& pw, int P) {
    int nsig = (int) pw.sig.size();
    res.vg.assign(nsig, std::vector<double>(P, NA_REAL));
    res.se_vg.assign(nsig, std::vector<double>(P, NA_REAL));
    res.h2.assign(nsig, std::vector<double>(P, NA_REAL));
    res.se_h2.assign(nsig, std::vector<double>(P, NA_REAL));
    res.fit.assign(P, NA_REAL);
    res.chr = pw.chr; res.w_start = pw.w_start;
    res.nL = pw.nL; res.nR = pw.nR; res.nC = pw.nC; res.mid_total = pw.mid_total;
    res.sig_cat = pw.sig_cat; res.sig_m = pw.sig_m;
}

// ---- HE (randomized MoM) with a genome-wide background component ----------
// Z  : SHARED probe block (n x B), same for every window.
// gw : precomputed genome-wide background (may be inactive).
static void he_part_compute(const PartWindow& pw, const Eigen::MatrixXd& Y,
                            const std::vector<double>& Vp, const GenoMat& Z,
                            const GwBackground& gw, bool se, PartResult& res) {
    int C = (int) pw.X.size(), n = (int) Y.rows(), P = (int) Y.cols();
    int nsig = (int) pw.sig.size();
    alloc_result(res, pw, P);
    if (C == 0) return;

    const bool useg = gw.active;
    const int gC  = C;                       // genome-wide component index (if useg)
    const int env = C + (useg ? 1 : 0);      // residual index
    const int dim = env + 1;
    const int B = (int) Z.cols();

    std::vector<double> M(C);
    for (int c = 0; c < C; ++c) M[c] = (double) pw.X[c].cols();

    GenoMat Yf = Y.cast<float>();
    std::vector<GenoMat> XtY(C);
    Eigen::MatrixXd T = Eigen::MatrixXd::Zero(dim, dim);
    for (int c = 0; c < C; ++c) {
        XtY[c] = pw.X[c].transpose() * Yf;
        double tr = pw.X[c].colwise().squaredNorm().cast<double>().sum() / M[c];   // tr(K_c) = n
        T(c, env) = tr; T(env, c) = tr;
    }
    T(env, env) = n;

    // K_c z_b for every window component (shared probes Z)
    std::vector<GenoMat> Gz(C);
    for (int c = 0; c < C; ++c) Gz[c] = (pw.X[c] * (pw.X[c].transpose() * Z)) / (float) M[c];
    for (int a = 0; a < C; ++a) for (int b = a; b < C; ++b) {
        double t = 0.0; for (int p = 0; p < B; ++p) t += (double) Gz[a].col(p).dot(Gz[b].col(p));
        T(a, b) = t / B; T(b, a) = T(a, b);
    }

    // genome-wide background: one extra row/col of T, using the precomputed U
    if (useg) {
        T(gC, env) = gw.trK; T(env, gC) = gw.trK;      // tr(K_gw) = n
        T(gC, gC)  = gw.trK2;                          // tr(K_gw^2)
        for (int a = 0; a < C; ++a) {                  // tr(K_a K_gw) ~ (1/B) sum_b K_a z_b . U_b
            double t = 0.0;
            for (int p = 0; p < B; ++p) t += (double) Gz[a].col(p).dot(gw.U.col(p));
            t /= B; T(a, gC) = t; T(gC, a) = t;
        }
    }

    Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> Tcod(T);

    // Reduced (window + residual) system for SE: the genome-wide background is
    // treated as a KNOWN offset (its sampling variance is negligible), so the
    // SE is computed exactly as in the no-gw estimator over C+1 components.
    Eigen::MatrixXd Tinv_red;
    if (se) {
        Eigen::MatrixXd Tred(C + 1, C + 1);
        Tred.topLeftCorner(C, C) = T.topLeftCorner(C, C);
        for (int c = 0; c < C; ++c) { Tred(c, C) = T(c, env); Tred(C, c) = T(env, c); }
        Tred(C, C) = n;
        Tinv_red = Tred.completeOrthogonalDecomposition().pseudoInverse();
    }

    for (int t = 0; t < P; ++t) {
        Eigen::VectorXd q(dim);
        for (int c = 0; c < C; ++c) q[c] = XtY[c].col(t).cast<double>().squaredNorm() / M[c];
        if (useg) q[gC] = gw.q[t];
        q[env] = Y.col(t).squaredNorm();
        Eigen::VectorXd sigma = Tcod.solve(q);
        res.fit[t] = sigma.dot(q);

        Eigen::MatrixXd Cov;
        if (se) {
            // sigmas for the reduced (window + residual) model: window comps
            // keep their solved values; residual index maps env -> C.
            Eigen::VectorXd sr(C + 1);
            for (int c = 0; c < C; ++c) sr[c] = sigma[c];
            sr[C] = sigma[env];

            GenoMat u = (float) sr[C] * Z;                 // Sigma z (window comps + I only)
            for (int c = 0; c < C; ++c) u.noalias() += (float)(sr[c] / M[c]) * (pw.X[c] * (pw.X[c].transpose() * Z));
            std::vector<GenoMat> skju(C + 1);
            for (int l = 0; l <= C; ++l) {
                GenoMat Klu;
                if (l == C) Klu = u; else Klu = (pw.X[l] * (pw.X[l].transpose() * u)) / (float) M[l];
                GenoMat sk = (float) sr[C] * Klu;
                for (int c = 0; c < C; ++c) sk.noalias() += (float)(sr[c] / M[c]) * (pw.X[c] * (pw.X[c].transpose() * Klu));
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
            Cov = Tinv_red * Covq * Tinv_red;
        }

        for (int k = 0; k < nsig; ++k) {
            int c = pw.sig[k];                             // window-component index (0..C-1)
            double vg = sigma[c];
            res.vg[k][t] = vg;
            res.h2[k][t] = (Vp[t] > 0) ? vg / Vp[t] : NA_REAL;
            if (se && Cov(c, c) > 0) {
                double sev = std::sqrt(Cov(c, c));
                res.se_vg[k][t] = sev;
                res.se_h2[k][t] = (Vp[t] > 0) ? sev / Vp[t] : NA_REAL;
            }
        }
    }
}

// --------------------------------------------------------------------------
// Progress bar (same throttled single-line style)
// --------------------------------------------------------------------------
static int count_total_windows(const RunContext& ctx, long W) {
    int total = 0;
    for (size_t ci = 0; ci < ctx.chr_order.size(); ++ci) {
        int c_lo = ctx.chr_lo[ci], c_hi = ctx.chr_hi[ci];
        long min_bp = ctx.wes_bim.bp[c_lo], max_bp = ctx.wes_bim.bp[c_hi];
        long first = (min_bp / W) * W;
        for (long w = first; w <= max_bp; w += W) {
            int iM0 = lower_index(ctx.wes_bim.bp, c_lo, c_hi, w);
            int iR0 = lower_index(ctx.wes_bim.bp, c_lo, c_hi, w + W);
            if (iR0 - iM0 > 0) ++total;
        }
    }
    return total;
}

static std::string fmt_time(double s) {
    if (s < 0 || !(s == s)) s = 0;
    int t = (int)(s + 0.5);
    char b[24];
    std::snprintf(b, sizeof(b), "%d:%02d", t / 60, t % 60);
    return std::string(b);
}

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
        if (!last && secs(tlast, now) < 0.2) return;
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

// --------------------------------------------------------------------------
// Driver (parallel over windows; streams a long-format table)
// --------------------------------------------------------------------------
struct PartParams {
    long W; int nmcmc; bool se;
    std::string out_file; int batch_size; int n_threads; unsigned seed;
    int gw_block;                          // genome-wide precompute SNP block
};

// Workers BUILD each window's components then run HE with the shared probes +
// genome-wide background.
struct HEWorker : public RcppParallel::Worker {
    const std::vector<RawWindow>& batch; const RunContext& ctx;
    const Eigen::MatrixXd& Y; const std::vector<double>& Vp;
    const GenoMat& Z; const GwBackground& gw;
    bool se; std::vector<PartResult>& out;
    HEWorker(const std::vector<RawWindow>& b, const RunContext& ctx, const Eigen::MatrixXd& Y,
             const std::vector<double>& Vp, const GenoMat& Z, const GwBackground& gw,
             bool se, std::vector<PartResult>& out)
        : batch(b), ctx(ctx), Y(Y), Vp(Vp), Z(Z), gw(gw), se(se), out(out) {}
    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t w = begin; w < end; ++w) {
            PartWindow pw; build_part_window(batch[w], ctx, pw);
            he_part_compute(pw, Y, Vp, Z, gw, se, out[w]);
        }
    }
};

static Rcpp::List part_driver(RunContext& ctx, const PartParams& pr) {
    const long W = pr.W; const int P = ctx.n_pheno; const Eigen::MatrixXd& Y = ctx.Y;
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

    // ONE shared probe block for the whole run (so the genome-wide U is
    // reusable in every window). Reproducible from `seed`.
    int B = pr.nmcmc;
    GenoMat Z(ctx.n_inds, B);
    { std::mt19937 rng(pr.seed); std::normal_distribution<float> nd(0.f, 1.f);
      for (int p = 0; p < B; ++p) for (int i = 0; i < ctx.n_inds; ++i) Z(i, p) = nd(rng); }

    // ONE streaming pass to precompute the genome-wide background (if common set).
    GwBackground gw = precompute_gw(ctx, Z, pr.gw_block);

    std::ofstream fout;
    bool tofile = !pr.out_file.empty();
    auto wr = [](std::ofstream& f, double v) { if (std::isnan(v)) f << "NA"; else f << v; };
    if (tofile) {
        fout.open(pr.out_file.c_str());
        if (!fout.is_open()) stop("Could not open out_file: " + pr.out_file);
        fout << "chr\tstart\tend\tn_left\tn_right\tn_common\tphenotype\tcategory\tm_c\tvg\tse_vg\th2\tse_h2\tenrichment\tfit\n";
    }

    int NC = ctx.n_cat + 1;
    std::vector< std::vector<double> > tot_vg(NC, std::vector<double>(P, 0.0));
    std::vector< std::vector<double> > tot_var(NC, std::vector<double>(P, 0.0));
    std::vector<double> tot_fit(P, 0.0);

    PartStream stream(ctx, W);
    long processed = 0; int n_win = 0;
    Progress prog; prog.start("HE-gw", count_total_windows(ctx, W));

    while (true) {
        std::vector<RawWindow> batch;
        RawWindow rw;
        while ((int) batch.size() < pr.batch_size && stream.next(rw)) batch.push_back(std::move(rw));
        if (batch.empty()) break;

        std::vector<PartResult> results(batch.size());
        { HEWorker wk(batch, ctx, Y, Vp, Z, gw, pr.se, results);
          RcppParallel::parallelFor(0, batch.size(), wk); }

        for (size_t b = 0; b < batch.size(); ++b) {
            const PartResult& r = results[b];
            int nsig = (int) r.sig_cat.size();
            for (int t = 0; t < P; ++t) {
                std::string ph = as<std::string>(ctx.trait_names[t]);
                if (!std::isnan(r.fit[t])) tot_fit[t] += r.fit[t];
                double mid_vg = 0.0;
                for (int k = 0; k < nsig; ++k) if (!std::isnan(r.vg[k][t])) mid_vg += r.vg[k][t];
                for (int k = 0; k < nsig; ++k) {
                    int cat = r.sig_cat[k]; int mc = r.sig_m[k];
                    double vg = r.vg[k][t];
                    double expected = (r.mid_total > 0) ? (double) mc / r.mid_total : NA_REAL;
                    double share = (mid_vg != 0.0 && !std::isnan(vg)) ? vg / mid_vg : NA_REAL;
                    double enr = (!std::isnan(share) && !std::isnan(expected) && expected > 0) ? share / expected : NA_REAL;
                    if (tofile) {
                        fout << r.chr << '\t' << r.w_start << '\t' << (r.w_start + W) << '\t'
                             << r.nL << '\t' << r.nR << '\t' << r.nC << '\t' << ph << '\t'
                             << ctx.cat_names[cat] << '\t' << mc << '\t';
                        wr(fout, vg); fout << '\t'; wr(fout, r.se_vg[k][t]); fout << '\t';
                        wr(fout, r.h2[k][t]); fout << '\t'; wr(fout, r.se_h2[k][t]); fout << '\t';
                        wr(fout, enr); fout << '\t'; wr(fout, r.fit[t]); fout << '\n';
                    }
                    if (!std::isnan(vg)) tot_vg[cat][t] += vg;
                    if (!std::isnan(r.se_vg[k][t])) tot_var[cat][t] += r.se_vg[k][t] * r.se_vg[k][t];
                }
            }
            ++n_win;
            prog.tick(r.chr, r.mid_total);
        }
        processed += (long) batch.size();
        if (tofile) fout.flush();
        Rcpp::checkUserInterrupt();
    }

    // genome-wide TOTAL rows (one per category per trait)
    if (tofile) {
        for (int t = 0; t < P; ++t) {
            std::string ph = as<std::string>(ctx.trait_names[t]);
            for (int c = 0; c < NC; ++c) {
                double vg = tot_vg[c][t], sev = std::sqrt(tot_var[c][t]);
                fout << "TOTAL\tNA\tNA\tNA\tNA\tNA\t" << ph << '\t' << ctx.cat_names[c] << "\tNA\t";
                wr(fout, vg); fout << '\t'; wr(fout, sev); fout << '\t';
                wr(fout, Vp[t] > 0 ? vg / Vp[t] : NA_REAL); fout << '\t';
                wr(fout, Vp[t] > 0 ? sev / Vp[t] : NA_REAL); fout << "\tNA";
                fout << '\t'; wr(fout, tot_fit[t]); fout << '\n';
            }
        }
        fout.close();
    }
    Rcout << "Windows estimated: " << n_win << "\n";

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

    Rcpp::NumericVector fit_out(P);
    for (int t = 0; t < P; ++t) fit_out[t] = tot_fit[t];
    fit_out.names() = ctx.trait_names;

    return List::create(_["genome_h2"] = gh2, _["genome_se"] = gse,
                        _["fit"] = fit_out,
                        _["categories"] = wrap(ctx.cat_names), _["trait_names"] = ctx.trait_names,
                        _["alpha"] = ctx.alpha, _["method"] = "HE",
                        _["gw_background"] = gw.active, _["gw_n_snps"] = (double) gw.M,
                        _["n_windows"] = n_win);
}

}  // end anonymous namespace

// ===========================================================================
// Export (external linkage; same arguments as he_sliding_window_part, plus an
// optional gw_block that only tunes the genome-wide precompute chunk size).
// ===========================================================================
// [[Rcpp::export]]
Rcpp::List he_sliding_window_part_gw(const std::string& filename,
                                     const SEXP pheno_mat,
                                     const IntegerMatrix snp_cat,
                                     const CharacterVector cat_names,
                                     double window_size = 1e6,
                                     double alpha = -1.0,
                                     Rcpp::Nullable<Rcpp::String> common_filename = R_NilValue,
                                     int nmcmc = 20,
                                     bool se = true,
                                     std::string out_file = "",
                                     int batch_size = 8,
                                     int n_threads = 0,
                                     int seed = 12345,
                                     int gw_block = 2048) {
    if (window_size <= 0) stop("window_size must be positive");
    RunContext ctx = setup_context(filename, pheno_mat, snp_cat, cat_names, alpha, common_filename);
    PartParams pr;
    pr.W = (long) window_size; pr.nmcmc = nmcmc; pr.se = se;
    pr.out_file = out_file; pr.batch_size = batch_size; pr.n_threads = n_threads;
    pr.seed = (unsigned) seed; pr.gw_block = gw_block;
    return part_driver(ctx, pr);
}
