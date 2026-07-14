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
#include <limits>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <cmath>

using namespace Rcpp;

// ===========================================================================
// Sliding-window local heritability, inspired by FLEX (Liu et al. 2025).
//
//   Tile the genome with a middle window of size W (default 1 Mb). For each
//   middle window we jointly fit THREE WES GRMs:
//        left  window  [w-W , w  )
//        middle window [w   , w+W)   <-- the component we record
//        right window  [w+W , w+2W)
//   The flanking components soak up "background h2" that leaks in through LD,
//   so the recorded middle component is the local signal conditional on its
//   neighbourhood (FLEX's "conditional" gene-vs-flank idea). Slide W right.
//
//   An optional COMMON-SNP GRM (separate PLINK fileset, restricted to the same
//   3W region) can be added as an extra background component, so WES variants
//   don't inflate h2 by tagging common variation.
//
// Two estimators over the SAME model V = sum_l s_l G_l + s_e I are provided:
//   * he_sliding_window()   - randomized Haseman-Elston / method-of-moments
//                             (this is what FLEX-h2 itself uses). Fast; also
//                             returns multi-trait genetic covariances.
//   * reml_sliding_window() - AI-REML (restricted maximum likelihood).
//                             Slower (dense O(N^3) per iteration) but the
//                             statistically efficient likelihood estimator the
//                             FLEX paper notes as the ideal-case reference.
//                             Univariate (per-trait h2).
//
// G_l = X_l X_l' / M_l  with X_l the column-standardized genotypes of window l.
// ===========================================================================

// --------------------------------------------------------------------------
// Shared helpers
// --------------------------------------------------------------------------

// Position info parsed from one PLINK .bim file.
struct BimInfo {
    std::vector<std::string> chr;   // chromosome label per SNP
    std::vector<long>        bp;    // base-pair position per SNP
    int n_snps;
};

// Parse chromosome (col 1) and bp (col 4) from a .bim file.
// Assumes SNPs are sorted by position within each chromosome (standard PLINK).
static BimInfo read_bim_positions(const std::string& bim_path) {
    BimInfo out;
    std::ifstream in(bim_path.c_str());
    if (!in.is_open()) stop("Could not open .bim file: " + bim_path);
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::istringstream ss(line);
        std::string chr, rsid, cm, bp;
        ss >> chr >> rsid >> cm >> bp;   // remaining a1/a2 ignored
        out.chr.push_back(chr);
        out.bp.push_back(std::stol(bp));
    }
    out.n_snps = (int) out.bp.size();
    return out;
}

// First index i in [lo, hi] with bp[i] >= value (contiguous, sorted region).
static int lower_index(const std::vector<long>& bp, int lo, int hi, long value) {
    std::vector<long>::const_iterator b = bp.begin();
    std::vector<long>::const_iterator it =
        std::lower_bound(b + lo, b + hi + 1, value);
    return (int)(it - b);
}

// Genotypes are stored in SINGLE precision (float) -- they are small integers
// (0/1/2, standardized to ~unit variance), and the big genotype matrix products
// dominate memory bandwidth at biobank N. float halves the bandwidth and
// doubles SIMD width (~2x). All *reductions* (means, variances, traces, dot
// products over N) are still accumulated in double to avoid loss of precision.
typedef Eigen::MatrixXf GenoMat;

// Column-standardize a genotype matrix in place: missing (== -1) -> column mean,
// then centre and scale to unit variance (denominator n-1); monomorphic columns
// -> 0. Mean/variance are accumulated in double; only the stored values are float.
static void standardize_columns(GenoMat& X) {
    const int n = (int) X.rows();
    for (int col = 0; col < X.cols(); ++col) {
        double sum = 0.0; int n_valid = 0;
        for (int i = 0; i < n; ++i) { float g = X(i, col); if (g != -1.f) { sum += g; ++n_valid; } }
        if (n_valid == 0) { X.col(col).setZero(); continue; }
        double mean = sum / n_valid;
        double sq_sum = 0.0;
        for (int i = 0; i < n; ++i) { double g = (X(i, col) == -1.f) ? mean : (double) X(i, col); sq_sum += (g - mean) * (g - mean); }
        double stddev = std::sqrt(sq_sum / (n - 1));
        if (stddev > 1e-10) {
            for (int i = 0; i < n; ++i) { double g = (X(i, col) == -1.f) ? mean : (double) X(i, col); X(i, col) = (float)((g - mean) / stddev); }
        } else {
            X.col(col).setZero();                    // monomorphic
        }
    }
}

// Read a contiguous SNP range from a .bed, subset to `keep_rows` individuals,
// return an (n_keep x n_cols) float matrix with missing preserved as -1.
static GenoMat read_geno_float(const std::string& bed_prefix,
                               int n_total_inds, int n_snps,
                               const std::vector<int>& keep_rows,
                               int snp_lo, int snp_hi) {
    int n_cols = snp_hi - snp_lo + 1;
    int n_keep = (int) keep_rows.size();
    GenoMat X(n_keep, n_cols);
    if (n_cols <= 0) return X;
    IntegerMatrix block = readBedBlock(
        bed_prefix + ".bed", n_total_inds, n_snps,
        0, n_total_inds - 1, snp_lo, snp_hi);
    for (int i = 0; i < n_keep; ++i)
        for (int j = 0; j < n_cols; ++j)
            X(i, j) = (float) block(keep_rows[i], j);
    return X;
}

// Standardized genotype components for one middle window.
struct WindowComponents {
    std::vector<std::string> label;      // "left","middle","right","common"
    std::vector<GenoMat>     X;          // standardized genotypes per component (float)
    int mid_comp;                          // index of the middle component (-1 = skip)
    int nL, nM, nR, nC;                    // SNP counts
};

// Read + standardize the genotypes of ONE window-sized cell: the SNPs of a
// fileset whose bp fall in [cell_start, cell_start + W) on chromosome index
// range [lo, hi]. Returns an (n_keep x m) matrix (m may be 0). This is the unit
// that the rolling WindowStream caches and shifts so each cell is read once.
static GenoMat read_cell(const std::string& prefix, int n_total, int n_snps,
                         const BimInfo& bim, const std::vector<int>& keep,
                         int lo, int hi, long cell_start, long W) {
    if (lo < 0 || hi < lo) return GenoMat(keep.size(), 0);
    int i0 = lower_index(bim.bp, lo, hi, cell_start);
    int i1 = lower_index(bim.bp, lo, hi, cell_start + W);   // one past
    int m = i1 - i0;
    if (m <= 0) return GenoMat(keep.size(), 0);
    GenoMat X = read_geno_float(prefix, n_total, n_snps, keep, i0, i1 - 1);
    standardize_columns(X);
    return X;
}

// Context bundled once and reused across all windows (avoids re-parsing files).
struct RunContext {
    std::string wes_prefix; int wes_n_total, wes_n_snps; BimInfo wes_bim;
    std::vector<int> geno_keep;
    Eigen::MatrixXd Y;                       // n_inds x n_pheno (WES order)
    CharacterVector trait_names;
    int n_inds, n_pheno;
    // chromosome index ranges (WES)
    std::vector<std::string> chr_order;
    std::vector<int> chr_lo, chr_hi;
    // common fileset
    bool use_common; std::string common_prefix;
    int common_n_total, common_n_snps; BimInfo common_bim;
    std::vector<int> common_keep;
};

// Common set-up (dimensions, id matching, phenotype ordering, chromosome
// ranges, optional common fileset) shared by both estimators.
static RunContext setup_context(const std::string& filename,
                                const SEXP pheno_mat,
                                Rcpp::Nullable<Rcpp::String> common_filename) {
    RunContext ctx;
    ctx.wes_prefix   = filename;
    ctx.wes_n_snps   = count_lines(filename + ".bim");
    List fam         = read_fam_file(filename);
    CharacterVector geno_iid = fam["iid"];
    ctx.wes_n_total  = geno_iid.size();
    Rcout << "WES individuals (.fam): " << ctx.wes_n_total << "\n";
    Rcout << "WES SNPs (.bim): "        << ctx.wes_n_snps  << "\n";

    ctx.wes_bim = read_bim_positions(filename + ".bim");
    if (ctx.wes_bim.n_snps != ctx.wes_n_snps)
        stop("Mismatch between .bim line count and parsed positions");

    // Phenotype
    Rcpp::NumericMatrix pheno;
    CharacterVector pheno_ids;
    if (Rf_isMatrix(pheno_mat) && !Rf_isNull(rownames(pheno_mat))) {
        pheno = as<NumericMatrix>(pheno_mat);
        pheno_ids = rownames(pheno_mat);
        SEXP cn = colnames(pheno_mat);
        if (!Rf_isNull(cn)) ctx.trait_names = cn;
    } else {
        stop("Phenotype must be a numeric matrix with IDs as rownames");
    }
    ctx.n_pheno = pheno.cols();
    if (ctx.trait_names.size() != ctx.n_pheno) {
        ctx.trait_names = CharacterVector(ctx.n_pheno);
        for (int j = 0; j < ctx.n_pheno; ++j)
            ctx.trait_names[j] = "trait" + std::to_string(j + 1);
    }

    IntegerVector match_idx = match(geno_iid, pheno_ids);
    std::vector<int> pheno_keep;
    for (int i = 0; i < match_idx.size(); ++i) {
        if (match_idx[i] != NA_INTEGER) {
            ctx.geno_keep.push_back(i);
            pheno_keep.push_back(match_idx[i] - 1);
        }
    }
    if (ctx.geno_keep.empty())
        stop("No overlapping individuals between genotype and phenotype");
    ctx.n_inds = (int) ctx.geno_keep.size();
    Rcout << "Individuals with complete data: " << ctx.n_inds << "\n";

    ctx.Y.resize(ctx.n_inds, ctx.n_pheno);
    for (int i = 0; i < ctx.n_inds; ++i)
        for (int j = 0; j < ctx.n_pheno; ++j)
            ctx.Y(i, j) = pheno(pheno_keep[i], j);

    CharacterVector analysis_iid(ctx.n_inds);
    for (int i = 0; i < ctx.n_inds; ++i)
        analysis_iid[i] = geno_iid[ctx.geno_keep[i]];

    // Chromosome index ranges (WES), preserving appearance order.
    for (int j = 0; j < ctx.wes_n_snps; ++j) {
        if (ctx.chr_order.empty() || ctx.wes_bim.chr[j] != ctx.chr_order.back()) {
            ctx.chr_order.push_back(ctx.wes_bim.chr[j]);
            ctx.chr_lo.push_back(j);
            ctx.chr_hi.push_back(j);
        } else {
            ctx.chr_hi.back() = j;
        }
    }

    // Optional common fileset
    ctx.use_common = common_filename.isNotNull();
    ctx.common_n_total = ctx.common_n_snps = 0;
    if (ctx.use_common) {
        ctx.common_prefix = as<std::string>(common_filename.get());
        ctx.common_n_snps = count_lines(ctx.common_prefix + ".bim");
        List cfam = read_fam_file(ctx.common_prefix);
        CharacterVector common_iid = cfam["iid"];
        ctx.common_n_total = common_iid.size();
        ctx.common_bim = read_bim_positions(ctx.common_prefix + ".bim");
        Rcout << "Common-SNP individuals (.fam): " << ctx.common_n_total << "\n";
        Rcout << "Common SNPs (.bim): "            << ctx.common_n_snps  << "\n";

        IntegerVector cmatch = match(analysis_iid, common_iid);
        ctx.common_keep.resize(ctx.n_inds);
        for (int i = 0; i < ctx.n_inds; ++i) {
            if (cmatch[i] == NA_INTEGER)
                stop("An analysis individual is missing from the common-SNP "
                     ".fam file; the common fileset must contain all analysis "
                     "individuals.");
            ctx.common_keep[i] = cmatch[i] - 1;
        }
    }
    return ctx;
}

// Common chromosome index range in the common fileset for a chromosome label.
static void common_chr_range(const RunContext& ctx, const std::string& chr,
                             int& cc_lo, int& cc_hi) {
    cc_lo = cc_hi = -1;
    if (!ctx.use_common) return;
    for (int j = 0; j < ctx.common_n_snps; ++j)
        if (ctx.common_bim.chr[j] == chr) {
            if (cc_lo == -1) cc_lo = j;
            cc_hi = j;
        }
}

// Rolling window iterator: yields the L/M/R (+ optional common) components for
// each middle window WITHOUT re-reading genotypes. It caches the three
// window-sized cells (g-1, g, g+1) and, when it advances one window to the
// right, shifts them (left <- middle, middle <- right) and reads only the new
// right cell. So each genotype cell is read from disk exactly once per pass,
// instead of ~3x (once as right, once as middle, once as left). The common-SNP
// cells are cached and shifted the same way. Shared by both HE and REML.
class WindowStream {
public:
    WindowStream(const RunContext& ctx, long W)
        : ctx_(ctx), W_(W), ci_(0), have_chrom_(false) {}

    // Fills wc/chr/w_start with the next middle window that has >=1 middle SNP.
    // Returns false when all chromosomes are exhausted.
    bool next(WindowComponents& wc, std::string& chr, long& w_start) {
        while (true) {
            if (!have_chrom_) {
                if (ci_ >= ctx_.chr_order.size()) return false;
                setup_chrom();
                have_chrom_ = true;
            }
            if (g_ > glast_) { ++ci_; have_chrom_ = false; continue; }

            bool produced = false;
            if (wesM_.cols() > 0) {                 // middle has SNPs -> a window
                assemble(wc);
                chr = ctx_.chr_order[ci_];
                w_start = g_ * W_;
                produced = true;
            }
            shift();                                 // advance to next middle cell
            if (produced) return true;
        }
    }

private:
    const RunContext& ctx_;
    long W_;
    size_t ci_;
    bool have_chrom_;
    int c_lo_, c_hi_, cc_lo_, cc_hi_;
    long g_, glast_;
    GenoMat wesL_, wesM_, wesR_;     // cells g-1, g, g+1 (WES)
    GenoMat comL_, comM_, comR_;     // cells g-1, g, g+1 (common)

    GenoMat wesCell(long cell) const {
        return read_cell(ctx_.wes_prefix, ctx_.wes_n_total, ctx_.wes_n_snps,
                         ctx_.wes_bim, ctx_.geno_keep, c_lo_, c_hi_, cell * W_, W_);
    }
    GenoMat comCell(long cell) const {
        if (!ctx_.use_common) return GenoMat(ctx_.n_inds, 0);
        return read_cell(ctx_.common_prefix, ctx_.common_n_total, ctx_.common_n_snps,
                         ctx_.common_bim, ctx_.common_keep, cc_lo_, cc_hi_, cell * W_, W_);
    }

    void setup_chrom() {
        const std::string& chr = ctx_.chr_order[ci_];
        c_lo_ = ctx_.chr_lo[ci_]; c_hi_ = ctx_.chr_hi[ci_];
        common_chr_range(ctx_, chr, cc_lo_, cc_hi_);
        long min_bp = ctx_.wes_bim.bp[c_lo_], max_bp = ctx_.wes_bim.bp[c_hi_];
        g_ = min_bp / W_;                    // first middle cell (grid aligned)
        glast_ = max_bp / W_;                // last cell that can be a middle
        wesL_ = wesCell(g_ - 1); wesM_ = wesCell(g_); wesR_ = wesCell(g_ + 1);
        comL_ = comCell(g_ - 1); comM_ = comCell(g_); comR_ = comCell(g_ + 1);
    }

    void shift() {
        ++g_;
        if (g_ > glast_) return;             // chromosome will roll over in next()
        wesL_ = std::move(wesM_); wesM_ = std::move(wesR_); wesR_ = wesCell(g_ + 1);
        comL_ = std::move(comM_); comM_ = std::move(comR_); comR_ = comCell(g_ + 1);
    }

    void assemble(WindowComponents& wc) {
        wc = WindowComponents();
        wc.mid_comp = -1;
        wc.nL = (int) wesL_.cols(); wc.nM = (int) wesM_.cols(); wc.nR = (int) wesR_.cols();
        if (wc.nL > 0) { wc.label.push_back("left");  wc.X.push_back(wesL_); }
        wc.label.push_back("middle"); wc.X.push_back(wesM_);
        wc.mid_comp = (int) wc.X.size() - 1;
        if (wc.nR > 0) { wc.label.push_back("right"); wc.X.push_back(wesR_); }

        int nC = (int) comL_.cols() + (int) comM_.cols() + (int) comR_.cols();
        wc.nC = nC;
        if (ctx_.use_common && nC > 0) {
            GenoMat Xc(ctx_.n_inds, nC);
            int off = 0;
            if (comL_.cols() > 0) { Xc.middleCols(off, comL_.cols()) = comL_; off += comL_.cols(); }
            if (comM_.cols() > 0) { Xc.middleCols(off, comM_.cols()) = comM_; off += comM_.cols(); }
            if (comR_.cols() > 0) { Xc.middleCols(off, comR_.cols()) = comR_; off += comR_.cols(); }
            wc.label.push_back("common");
            wc.X.push_back(Xc);
        }
    }
};

// Build an empty windows data.frame skeleton and finalize it with row.names.
static void finalize_dataframe(List& df, int n_win) {
    df.attr("class") = "data.frame";
    IntegerVector rn(n_win);
    for (int i = 0; i < n_win; ++i) rn[i] = i + 1;
    df.attr("row.names") = rn;
}

// Count windows that will actually be processed (middle window has >=1 SNP),
// without any genotype I/O, so we can show a percentage / ETA.
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
        Rcout << "[" << tag << "] scanning " << total
              << " windows...\n";
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

// Parameters passed to the shared parallel + streaming driver (defined at the
// end of the file, after the per-window compute functions).
struct DriverParams {
    long W; int method;                    // method: 0 = HE, 1 = REML
    int nmcmc; int max_iter; double tol;
    bool se;                               // compute standard errors?
    std::string out_file; int batch_size; int n_threads; unsigned seed;
};
static Rcpp::List sliding_window_driver(RunContext& ctx, const DriverParams& pr);

// ===========================================================================
// Estimator 1: randomized Haseman-Elston (method of moments). Parallelized over
// windows; streams per-window h2 + SE to `out_file` while running.
// ===========================================================================
// [[Rcpp::export]]
Rcpp::List he_sliding_window(const std::string& filename,
                             const SEXP pheno_mat,
                             double window_size = 1e6,
                             Rcpp::Nullable<Rcpp::String> common_filename = R_NilValue,
                             int nmcmc = 20,
                             bool se = true,
                             std::string out_file = "",
                             int batch_size = 64,
                             int n_threads = 0,
                             int seed = 12345) {
    if (window_size <= 0) stop("window_size must be positive");
    RunContext ctx = setup_context(filename, pheno_mat, common_filename);
    DriverParams pr;
    pr.W = (long) window_size; pr.method = 0;
    pr.nmcmc = nmcmc; pr.max_iter = 0; pr.tol = 0.0; pr.se = se;
    pr.out_file = out_file; pr.batch_size = batch_size;
    pr.n_threads = n_threads; pr.seed = (unsigned) seed;
    return sliding_window_driver(ctx, pr);
}

// ===========================================================================
// Estimator 2: AI-REML (restricted maximum likelihood), univariate per trait.
//
// Scalability note: a naive REML forms the N x N covariance V = sum_l s_l G_l
// + s_e I and inverts it (O(N^3), O(N^2) memory) every iteration -- impossible
// at biobank N (an N x N matrix at N=200k is ~320 GB). Because every GRM here
// is low rank, G_l = X_l X_l' / M_l, we use the Woodbury identity to run the
// ENTIRE AI-REML in the SNP dimension K = sum_l M_l (total SNPs in the 3-window
// region) instead of N. All heavy linear algebra is on K x K matrices; the
// only N-sized work is the one-off Gram matrix Z'Z (O(N K^2)) and a few
// matrix-vector products (O(N K)). No N x N matrix is ever formed.
//
// Let Z = [X_1 | ... | X_C] (N x K), w_l = s_l / M_l, Ws = diag(sqrt(w)) over
// columns, A = Z Ws so V = s_e I + A A'. With S = I_K + A'A / s_e:
//   V^-1        = (1/s_e) I - (1/s_e^2) A S^-1 A'
//   log|V|      = N log(s_e) + log|S|
// and every REML quantity (log-likelihood, scores, Average-Information matrix)
// reduces to K-space expressions in the precomputed Gram = Z'Z. Verified to
// match the dense N-space AI-REML to ~1e-14 on logL, score, AI and to ~1e-15
// on the fitted variance components.
// ===========================================================================

// Precomputed per-window sufficient statistics shared across traits/iterations.
struct RemlWindow {
    Eigen::MatrixXd Z;        // N x K standardized genotypes [X_1 | ... | X_C]
    Eigen::MatrixXd Gram;     // K x K  = Z'Z
    Eigen::VectorXd Zt1;      // K      = Z' 1  (column sums)
    std::vector<int> off, m;  // per-component column offset and size (M_l)
    int C, K;
};

struct RemlEval {
    double logL;
    Eigen::VectorXd score;    // length C+1 (residual last)
    Eigen::MatrixXd AI;       // (C+1) x (C+1)
    bool ok;
};

// Column scale vector Ws (sqrt of s_l / M_l per component block).
static Eigen::VectorXd reml_wsqrt(const RemlWindow& W, const Eigen::VectorXd& s) {
    Eigen::VectorXd Ws(W.K);
    for (int l = 0; l < W.C; ++l) {
        double w = s[l] / W.m[l];
        double sq = std::sqrt(w > 0 ? w : 0.0);
        for (int k = 0; k < W.m[l]; ++k) Ws[W.off[l] + k] = sq;
    }
    return Ws;
}

// Cheap path: REML log-likelihood only (used inside step-halving).
static bool reml_logL(const RemlWindow& W, const Eigen::VectorXd& Zty,
                      const Eigen::VectorXd& y, const Eigen::VectorXd& s,
                      double& logL) {
    double se = s[W.C];
    if (se <= 0) return false;
    Eigen::VectorXd Ws = reml_wsqrt(W, s);
    Eigen::MatrixXd S = (Ws.asDiagonal() * W.Gram * Ws.asDiagonal()) / se;
    S.diagonal().array() += 1.0;
    Eigen::LLT<Eigen::MatrixXd> llt(S);
    if (llt.info() != Eigen::Success) return false;
    double logdetS = 2.0 * llt.matrixLLT().diagonal().array().log().sum();

    Eigen::VectorXd g1 = Ws.cwiseProduct(W.Zt1), gy = Ws.cwiseProduct(Zty);
    Eigen::VectorXd v1 = Ws.cwiseProduct(llt.solve(g1));
    Eigen::VectorXd vy = Ws.cwiseProduct(llt.solve(gy));
    Eigen::VectorXd Vi1 = Eigen::VectorXd::Ones(W.Z.rows()) / se - (W.Z * v1) / (se * se);
    Eigen::VectorXd Viy = y / se - (W.Z * vy) / (se * se);
    double c = Vi1.sum();
    Eigen::VectorXd r = Viy - Vi1 * (Viy.sum() / c);     // r = P y
    double logdetV = (double) W.Z.rows() * std::log(se) + logdetS;
    logL = -0.5 * (logdetV + std::log(c) + y.dot(r));
    return true;
}

// Full path: log-likelihood, score vector, and Average-Information matrix.
static RemlEval reml_eval(const RemlWindow& W, const Eigen::VectorXd& Zty,
                          const Eigen::VectorXd& y, const Eigen::VectorXd& s) {
    RemlEval ev; ev.ok = false;
    int C = W.C, K = W.K, p = C + 1;
    double se = s[C];
    if (se <= 0) return ev;
    double Nn = (double) W.Z.rows();

    Eigen::VectorXd Ws = reml_wsqrt(W, s);
    Eigen::MatrixXd S = (Ws.asDiagonal() * W.Gram * Ws.asDiagonal()) / se;
    S.diagonal().array() += 1.0;
    Eigen::LLT<Eigen::MatrixXd> llt(S);
    if (llt.info() != Eigen::Success) return ev;
    Eigen::MatrixXd Sinv = llt.solve(Eigen::MatrixXd::Identity(K, K));
    double logdetS = 2.0 * llt.matrixLLT().diagonal().array().log().sum();

    Eigen::MatrixXd Gw  = W.Gram * Ws.asDiagonal();       // Gw_{jk} = G_{jk} Ws_k
    Eigen::MatrixXd R   = Sinv * Gw.transpose();
    Eigen::MatrixXd ZVZ = W.Gram / se - (Gw * R) / (se * se);   // Z' V^-1 Z

    Eigen::VectorXd g1 = Ws.cwiseProduct(W.Zt1), gy = Ws.cwiseProduct(Zty);
    Eigen::VectorXd s1 = Sinv * g1, sy = Sinv * gy;
    Eigen::VectorXd ZV1 = W.Zt1 / se - (Gw * s1) / (se * se);   // Z' V^-1 1

    Eigen::VectorXd Vi1 = Eigen::VectorXd::Ones((int)Nn) / se
                          - (W.Z * Ws.cwiseProduct(s1)) / (se * se);
    Eigen::VectorXd Viy = y / se - (W.Z * Ws.cwiseProduct(sy)) / (se * se);
    double c = Vi1.sum();
    Eigen::VectorXd r = Viy - Vi1 * (Viy.sum() / c);           // r = P y
    Eigen::VectorXd Zr = W.Z.transpose() * r;

    // Pr = P (P y)  (needed for AI entries involving the residual component)
    Eigen::VectorXd a_r = Ws.cwiseProduct(Zr);
    Eigen::VectorXd Vir = r / se - (W.Z * Ws.cwiseProduct(Sinv * a_r)) / (se * se);
    Eigen::VectorXd Pr  = Vir - Vi1 * (Vir.sum() / c);
    Eigen::VectorXd ZPr = W.Z.transpose() * Pr;

    double trVinv = Nn / se - (K - Sinv.trace()) / se;
    double trP    = trVinv - Vi1.dot(Vi1) / c;
    Eigen::MatrixXd Q = ZVZ - (ZV1 * ZV1.transpose()) / c;      // Z' P Z

    Eigen::VectorXd score(p);
    for (int l = 0; l < C; ++l) {
        double trVG = 0, xv1 = 0, pl2 = 0;
        for (int k = 0; k < W.m[l]; ++k) {
            trVG += ZVZ(W.off[l] + k, W.off[l] + k);
            double z1 = ZV1[W.off[l] + k]; xv1 += z1 * z1;
            double zr = Zr[W.off[l] + k];  pl2 += zr * zr;
        }
        trVG /= W.m[l]; xv1 /= W.m[l]; pl2 /= W.m[l];
        double trPG = trVG - xv1 / c;                          // tr(P G_l)
        score[l] = -0.5 * (trPG - pl2);                        // pl2 = y'P G_l P y
    }
    score[C] = -0.5 * (trP - r.dot(r));

    Eigen::MatrixXd AI(p, p);
    for (int a = 0; a < C; ++a)
        for (int b = a; b < C; ++b) {
            Eigen::VectorXd pa = Zr.segment(W.off[a], W.m[a]);
            Eigen::VectorXd pb = Zr.segment(W.off[b], W.m[b]);
            Eigen::MatrixXd Qs = Q.block(W.off[a], W.off[b], W.m[a], W.m[b]);
            double v = 0.5 / ((double) W.m[a] * W.m[b]) * pa.dot(Qs * pb);
            AI(a, b) = v; AI(b, a) = v;
        }
    for (int l = 0; l < C; ++l) {
        Eigen::VectorXd pl = Zr.segment(W.off[l], W.m[l]);
        Eigen::VectorXd zp = ZPr.segment(W.off[l], W.m[l]);
        double v = 0.5 / W.m[l] * pl.dot(zp);
        AI(l, C) = v; AI(C, l) = v;
    }
    AI(C, C) = 0.5 * r.dot(Pr);

    double logdetV = Nn * std::log(se) + logdetS;
    ev.logL = -0.5 * (logdetV + std::log(c) + y.dot(r));
    ev.score = score; ev.AI = AI; ev.ok = true;
    return ev;
}

// AI-REML for one trait using the low-rank window statistics. Returns variance
// components `s` (length C+1, residual last), convergence flag and iterations.
// Negative components are floored to keep S positive definite.
static void reml_fit(const RemlWindow& W, const Eigen::VectorXd& Zty,
                     const Eigen::VectorXd& y, int max_iter, double tol,
                     Eigen::VectorXd& s, bool& converged, int& iters,
                     Eigen::MatrixXd& AI_out) {
    int p = W.C + 1;
    double vy = (y.array() - y.mean()).square().sum() / (y.size() - 1);
    double floor_v = 1e-6 * vy;

    s = Eigen::VectorXd::Constant(p, vy / p);       // start: split variance
    converged = false; iters = 0;
    AI_out = Eigen::MatrixXd::Zero(p, p);

    RemlEval cur = reml_eval(W, Zty, y, s);
    if (!cur.ok) {                                  // nudge if start not PD
        s.setConstant(vy / p); s[W.C] = vy;
        cur = reml_eval(W, Zty, y, s);
        if (!cur.ok) return;
    }
    AI_out = cur.AI;

    for (int it = 1; it <= max_iter; ++it) {
        iters = it;
        Eigen::VectorXd delta = cur.AI.ldlt().solve(cur.score);

        // Step control via the cheap log-likelihood path.
        double step = 1.0; bool accepted = false; Eigen::VectorXd s_acc;
        for (int h = 0; h < 25; ++h) {
            Eigen::VectorXd s_try = s + step * delta;
            for (int l = 0; l < p; ++l) if (s_try[l] < floor_v) s_try[l] = floor_v;
            double ll;
            if (reml_logL(W, Zty, y, s_try, ll) && ll >= cur.logL - 1e-8) {
                s_acc = s_try; accepted = true; break;
            }
            step *= 0.5;
        }
        if (!accepted) break;

        s = s_acc;
        RemlEval full = reml_eval(W, Zty, y, s);
        if (!full.ok) break;
        AI_out = full.AI;
        double dlogL = full.logL - cur.logL;
        cur = full;
        if (std::abs(dlogL) < tol) { converged = true; break; }
    }
}

// ===========================================================================
// Estimator 2: AI-REML. Parallelized over windows; streams h2 + SE to file.
// ===========================================================================
// [[Rcpp::export]]
Rcpp::List reml_sliding_window(const std::string& filename,
                               const SEXP pheno_mat,
                               double window_size = 1e6,
                               Rcpp::Nullable<Rcpp::String> common_filename = R_NilValue,
                               int max_iter = 100,
                               double tol = 1e-4,
                               bool se = true,
                               std::string out_file = "",
                               int batch_size = 16,
                               int n_threads = 0,
                               int seed = 12345) {
    if (window_size <= 0) stop("window_size must be positive");
    RunContext ctx = setup_context(filename, pheno_mat, common_filename);
    DriverParams pr;
    pr.W = (long) window_size; pr.method = 1;
    pr.nmcmc = 0; pr.max_iter = max_iter; pr.tol = tol; pr.se = se;
    pr.out_file = out_file; pr.batch_size = batch_size;
    pr.n_threads = n_threads; pr.seed = (unsigned) seed;
    return sliding_window_driver(ctx, pr);
}

// ===========================================================================
// Per-window compute kernels (thread-safe: pure Eigen/STL, no R API) + the
// RcppParallel worker and the batched, streaming driver.
// ===========================================================================

// One window's results (all traits), filled by a worker thread.
struct WindowResult {
    int nL, nM, nR, nC;
    std::vector<double> vg, h2, se_vg, se_h2;
    std::vector<int> conv, iters;            // REML only
};

// Apply Sigma = sum_c (sig[c]/M[c]) X_c X_c' + sig[C] I to a block of probe
// vectors (columns), in float. All big products are float GEMMs.
static inline GenoMat apply_Sigma_mat(const std::vector<GenoMat>& X,
                                      const std::vector<double>& M,
                                      const Eigen::VectorXd& sig, const GenoMat& Zb) {
    int C = (int) X.size();
    GenoMat o = (float) sig[C] * Zb;
    for (int c = 0; c < C; ++c)
        o.noalias() += (float)(sig[c] / M[c]) * (X[c] * (X[c].transpose() * Zb));
    return o;
}
// Apply K_c = X_c X_c'/M_c (identity for c == C) to a block of probe vectors.
static inline GenoMat apply_K_mat(const std::vector<GenoMat>& X,
                                  const std::vector<double>& M,
                                  int c, const GenoMat& Zb) {
    if (c == (int) X.size()) return Zb;
    return (X[c] * (X[c].transpose() * Zb)) / (float) M[c];
}

// ---- HE (randomized method of moments) with moment-based standard errors ----
// Thread-safe: probes use a per-window std::mt19937 seeded from `seed`. All
// genotype products are single-precision GEMMs with the B probes batched into
// one matrix; the moment matrix, q, and Cov(q) are accumulated in double.
static void he_window_compute(const WindowComponents& wc, const Eigen::MatrixXd& Y,
                              int nmcmc, unsigned seed, bool se, WindowResult& res) {
    int C = (int) wc.X.size(), env = C, n = (int) Y.rows(), P = (int) Y.cols();
    int mid = wc.mid_comp;
    res.nL = wc.nL; res.nM = wc.nM; res.nR = wc.nR; res.nC = wc.nC;
    res.vg.assign(P, NA_REAL); res.h2.assign(P, NA_REAL);
    res.se_vg.assign(P, NA_REAL); res.se_h2.assign(P, NA_REAL);

    std::vector<double> M(C);
    for (int c = 0; c < C; ++c) M[c] = (double) wc.X[c].cols();

    GenoMat Yf = Y.cast<float>();                        // phenotype in float
    std::vector<GenoMat> XtY(C);
    Eigen::MatrixXd T = Eigen::MatrixXd::Zero(C + 1, C + 1);
    for (int c = 0; c < C; ++c) {
        XtY[c] = wc.X[c].transpose() * Yf;               // m_c x P (float)
        double tr = wc.X[c].colwise().squaredNorm().cast<double>().sum() / M[c];  // tr(G_c) exact, double
        T(c, env) = tr; T(env, c) = tr;
    }
    T(env, env) = n;

    // Batched probes: one N x B float matrix -> two GEMMs per component.
    int B = nmcmc;
    GenoMat Z(n, B);
    std::mt19937 rng(seed);
    std::normal_distribution<float> nd(0.f, 1.f);
    for (int p = 0; p < B; ++p) for (int i = 0; i < n; ++i) Z(i, p) = nd(rng);
    std::vector<GenoMat> Gz(C);
    for (int c = 0; c < C; ++c) Gz[c] = (wc.X[c] * (wc.X[c].transpose() * Z)) / (float) M[c];
    for (int a = 0; a < C; ++a)
        for (int b = a; b < C; ++b) {
            double t = 0.0;
            for (int p = 0; p < B; ++p) t += (double) Gz[a].col(p).dot(Gz[b].col(p));
            T(a, b) = t / B; T(b, a) = T(a, b);          // tr(K_a K_b), double-accumulated
        }
    Eigen::MatrixXd Tinv = T.inverse();

    for (int t = 0; t < P; ++t) {
        Eigen::VectorXd q(C + 1);
        for (int c = 0; c < C; ++c) q[c] = XtY[c].col(t).cast<double>().squaredNorm() / M[c];
        q[env] = Y.col(t).squaredNorm();
        Eigen::VectorXd sigma = Tinv * q;
        double tot = sigma.sum(), vg = sigma[mid];
        res.vg[t] = vg;
        res.h2[t] = (std::abs(tot) > 1e-12) ? vg / tot : NA_REAL;
        if (!se) continue;                               // skip the expensive SE pass

        // SE: Cov(sigma) = Tinv Cov(q) Tinv, Cov(q)_ij = 2 tr(K_i Sigma K_j Sigma).
        // Reuse the probe block Z and Gz (= K_i z). Batched float; double sums.
        GenoMat u = apply_Sigma_mat(wc.X, M, sigma, Z);          // Sigma z
        std::vector<GenoMat> skju(C + 1);
        for (int l = 0; l <= C; ++l)
            skju[l] = apply_Sigma_mat(wc.X, M, sigma, apply_K_mat(wc.X, M, l, u));
        Eigen::MatrixXd Covq(C + 1, C + 1);
        for (int i = 0; i <= C; ++i) {
            const GenoMat& kizi = (i < C) ? Gz[i] : Z;          // K_i z (identity for residual)
            for (int j = 0; j <= C; ++j) {
                double v = 0.0;
                for (int p = 0; p < B; ++p) v += (double) kizi.col(p).dot(skju[j].col(p));
                Covq(i, j) = v;
            }
        }
        Covq = (Covq / B) * 2.0;
        Covq = (0.5 * (Covq + Covq.transpose())).eval();
        Eigen::MatrixXd Cov = Tinv * Covq * Tinv;
        if (Cov(mid, mid) > 0) res.se_vg[t] = std::sqrt(Cov(mid, mid));
        if (std::abs(tot) > 1e-12) {
            double h2 = vg / tot;
            Eigen::VectorXd g(C + 1);
            for (int k = 0; k <= C; ++k) g[k] = ((k == mid ? 1.0 : 0.0) - h2) / tot;
            double v = g.dot(Cov * g);
            if (v > 0) res.se_h2[t] = std::sqrt(v);
        }
    }
}

// ---- AI-REML with SE from the Average-Information matrix (delta method) ----
static void reml_window_compute(const WindowComponents& wc, const Eigen::MatrixXd& Y,
                                int max_iter, double tol, bool se, WindowResult& res) {
    int C = (int) wc.X.size(), n = (int) Y.rows(), P = (int) Y.cols();
    int mid = wc.mid_comp;
    res.nL = wc.nL; res.nM = wc.nM; res.nR = wc.nR; res.nC = wc.nC;
    res.vg.assign(P, NA_REAL); res.h2.assign(P, NA_REAL);
    res.se_vg.assign(P, NA_REAL); res.se_h2.assign(P, NA_REAL);
    res.conv.assign(P, 0); res.iters.assign(P, 0);

    RemlWindow rw;
    rw.C = C; rw.off.resize(C); rw.m.resize(C);
    int K = 0;
    for (int c = 0; c < C; ++c) { rw.off[c] = K; rw.m[c] = (int) wc.X[c].cols(); K += rw.m[c]; }
    rw.K = K;

    // Gram Z'Z (the O(N K^2) bottleneck): concatenate float components and use a
    // symmetric rank update (half the flops) in single precision, then cast the
    // small K x K result to double. The solver runs in double on a double Z.
    GenoMat Zf(n, K);
    for (int c = 0; c < C; ++c) Zf.middleCols(rw.off[c], rw.m[c]) = wc.X[c];
    GenoMat Gf = GenoMat::Zero(K, K);
    Gf.selfadjointView<Eigen::Upper>().rankUpdate(Zf.transpose());   // Gf = Zf' Zf (upper)
    rw.Gram = GenoMat(Gf.selfadjointView<Eigen::Upper>()).cast<double>();
    rw.Z    = Zf.cast<double>();
    rw.Zt1  = rw.Z.colwise().sum().transpose();

    for (int t = 0; t < P; ++t) {
        Eigen::VectorXd y = Y.col(t);
        Eigen::VectorXd Zty = rw.Z.transpose() * y;
        Eigen::VectorXd s; bool conv; int iters; Eigen::MatrixXd AI;
        reml_fit(rw, Zty, y, max_iter, tol, s, conv, iters, AI);
        double tot = s.sum(), vg = s[mid];
        res.vg[t] = vg; res.conv[t] = conv ? 1 : 0; res.iters[t] = iters;
        res.h2[t] = (conv && std::abs(tot) > 1e-12) ? vg / tot : NA_REAL;
        if (se && conv) {
            Eigen::MatrixXd Cov = AI.inverse();          // asymptotic Cov(sigma)
            if (Cov(mid, mid) > 0) res.se_vg[t] = std::sqrt(Cov(mid, mid));
            if (std::abs(tot) > 1e-12) {
                double h2 = vg / tot; int p = (int) s.size();
                Eigen::VectorXd g(p);
                for (int k = 0; k < p; ++k) g[k] = ((k == mid ? 1.0 : 0.0) - h2) / tot;
                double v = g.dot(Cov * g);
                if (v > 0) res.se_h2[t] = std::sqrt(v);
            }
        }
    }
}

// RcppParallel worker: computes a contiguous chunk of the current batch.
struct WindowWorker : public RcppParallel::Worker {
    const std::vector<WindowComponents>& batch;
    const Eigen::MatrixXd& Y;
    const DriverParams& pr;
    unsigned seed_offset;
    std::vector<WindowResult>& out;
    WindowWorker(const std::vector<WindowComponents>& batch, const Eigen::MatrixXd& Y,
                 const DriverParams& pr, unsigned seed_offset, std::vector<WindowResult>& out)
        : batch(batch), Y(Y), pr(pr), seed_offset(seed_offset), out(out) {}
    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t w = begin; w < end; ++w) {
            if (pr.method == 0)
                he_window_compute(batch[w], Y, pr.nmcmc,
                                  pr.seed + seed_offset + (unsigned) w, pr.se, out[w]);
            else
                reml_window_compute(batch[w], Y, pr.max_iter, pr.tol, pr.se, out[w]);
        }
    }
};

// Batched, streaming, parallel driver shared by both estimators.
static Rcpp::List sliding_window_driver(RunContext& ctx, const DriverParams& pr) {
    const long W = pr.W;
    const int P = ctx.n_pheno;
    const Eigen::MatrixXd& Y = ctx.Y;

    Eigen::setNbThreads(1);                              // RcppParallel owns threads
    // setThreadOptions() is an R-level function, not part of the C++ API. The
    // RcppParallel backend reads the thread count from this env var instead.
    if (pr.n_threads > 0) {
        std::string nt = std::to_string(pr.n_threads);
#ifdef _WIN32
        _putenv_s("RCPP_PARALLEL_NUM_THREADS", nt.c_str());
#else
        setenv("RCPP_PARALLEL_NUM_THREADS", nt.c_str(), 1);
#endif
    }

    std::vector<double> Vp(P);                           // phenotypic variance / trait
    for (int t = 0; t < P; ++t) {
        double m = Y.col(t).mean();
        Vp[t] = (Y.col(t).array() - m).square().sum() / (Y.rows() - 1);
    }

    const bool tofile = !pr.out_file.empty();
    std::ofstream fout;
    auto wr = [](std::ofstream& f, double v) { if (std::isnan(v)) f << "NA"; else f << v; };
    if (tofile) {
        fout.open(pr.out_file.c_str());
        if (!fout.is_open()) stop("Could not open out_file for writing: " + pr.out_file);
        fout << "chr\tstart\tend\tn_left\tn_middle\tn_right";
        if (ctx.use_common) fout << "\tn_common";
        for (int t = 0; t < P; ++t) {
            std::string nm = as<std::string>(ctx.trait_names[t]);
            fout << "\tvg_" << nm << "\tse_vg_" << nm << "\th2_" << nm << "\tse_h2_" << nm;
            if (pr.method == 1) fout << "\tconv_" << nm << "\tn_iter_" << nm;
        }
        fout << "\n";
    }

    std::vector<std::string> a_chr;
    std::vector<double> a_start, a_end;
    std::vector<int> a_nL, a_nM, a_nR, a_nC;
    std::vector< std::vector<double> > a_vg(P), a_h2(P), a_sevg(P), a_seh2(P);
    std::vector< std::vector<int> > a_conv(P), a_iter(P);
    std::vector<double> tot_vg(P, 0.0), tot_var(P, 0.0);   // sum vg, sum se_vg^2

    Progress prog; prog.start(pr.method == 0 ? "HE" : "REML", count_total_windows(ctx, W));
    WindowStream stream(ctx, W);
    long processed = 0;

    while (true) {
        std::vector<WindowComponents> comps;
        std::vector<std::string> bchr; std::vector<long> bws;
        WindowComponents wc; std::string chr; long ws;
        while ((int) comps.size() < pr.batch_size && stream.next(wc, chr, ws)) {
            comps.push_back(std::move(wc)); bchr.push_back(chr); bws.push_back(ws);
        }
        if (comps.empty()) break;

        std::vector<WindowResult> results(comps.size());
        WindowWorker worker(comps, Y, pr, (unsigned) processed, results);
        RcppParallel::parallelFor(0, comps.size(), worker);

        for (size_t k = 0; k < comps.size(); ++k) {
            const WindowResult& r = results[k];
            prog.tick(bchr[k], r.nM);
            if (tofile) {
                fout << bchr[k] << '\t' << bws[k] << '\t' << (bws[k] + W)
                     << '\t' << r.nL << '\t' << r.nM << '\t' << r.nR;
                if (ctx.use_common) fout << '\t' << r.nC;
                for (int t = 0; t < P; ++t) {
                    fout << '\t'; wr(fout, r.vg[t]);
                    fout << '\t'; wr(fout, r.se_vg[t]);
                    fout << '\t'; wr(fout, r.h2[t]);
                    fout << '\t'; wr(fout, r.se_h2[t]);
                    if (pr.method == 1) fout << '\t' << r.conv[t] << '\t' << r.iters[t];
                }
                fout << '\n';
            }
            a_chr.push_back(bchr[k]); a_start.push_back((double) bws[k]);
            a_end.push_back((double)(bws[k] + W));
            a_nL.push_back(r.nL); a_nM.push_back(r.nM); a_nR.push_back(r.nR); a_nC.push_back(r.nC);
            for (int t = 0; t < P; ++t) {
                a_vg[t].push_back(r.vg[t]); a_h2[t].push_back(r.h2[t]);
                a_sevg[t].push_back(r.se_vg[t]); a_seh2[t].push_back(r.se_h2[t]);
                if (pr.method == 1) { a_conv[t].push_back(r.conv[t]); a_iter[t].push_back(r.iters[t]); }
                if (!std::isnan(r.vg[t]))    tot_vg[t]  += r.vg[t];
                if (!std::isnan(r.se_vg[t])) tot_var[t] += r.se_vg[t] * r.se_vg[t];
            }
        }
        processed += (long) comps.size();
        if (tofile) fout.flush();
        Rcpp::checkUserInterrupt();
    }

    // Total = sum of middle-window genetic variance over the genome; SE from
    // summed per-window variances (windows treated as ~independent segments).
    std::vector<double> tot_h2(P), tot_h2_se(P), tot_se_vg(P);
    for (int t = 0; t < P; ++t) {
        tot_se_vg[t] = std::sqrt(tot_var[t]);
        tot_h2[t]    = (Vp[t] > 0) ? tot_vg[t]    / Vp[t] : NA_REAL;
        tot_h2_se[t] = (Vp[t] > 0) ? tot_se_vg[t] / Vp[t] : NA_REAL;
    }
    if (tofile) {
        fout << "TOTAL\tNA\tNA\tNA\tNA\tNA";
        if (ctx.use_common) fout << "\tNA";
        for (int t = 0; t < P; ++t) {
            fout << '\t'; wr(fout, tot_vg[t]);
            fout << '\t'; wr(fout, tot_se_vg[t]);
            fout << '\t'; wr(fout, tot_h2[t]);
            fout << '\t'; wr(fout, tot_h2_se[t]);
            if (pr.method == 1) fout << "\tNA\tNA";
        }
        fout << "\n";
        fout.close();
    }

    int n_win = (int) a_chr.size();
    Rcout << "Windows estimated: " << n_win << "\n";

    List df;
    df["chr"] = wrap(a_chr); df["start"] = wrap(a_start); df["end"] = wrap(a_end);
    df["n_left"] = wrap(a_nL); df["n_middle"] = wrap(a_nM); df["n_right"] = wrap(a_nR);
    if (ctx.use_common) df["n_common"] = wrap(a_nC);
    for (int t = 0; t < P; ++t) {
        std::string nm = as<std::string>(ctx.trait_names[t]);
        df["vg_" + nm]    = wrap(a_vg[t]);
        df["se_vg_" + nm] = wrap(a_sevg[t]);
        df["h2_" + nm]    = wrap(a_h2[t]);
        df["se_h2_" + nm] = wrap(a_seh2[t]);
        if (pr.method == 1) {
            df["converged_" + nm] = wrap(a_conv[t]);
            df["n_iter_" + nm]    = wrap(a_iter[t]);
        }
    }
    finalize_dataframe(df, n_win);

    List tot;
    tot["trait"]        = ctx.trait_names;
    tot["total_vg"]     = wrap(tot_vg);
    tot["total_se_vg"]  = wrap(tot_se_vg);
    tot["total_h2"]     = wrap(tot_h2);
    tot["total_h2_se"]  = wrap(tot_h2_se);
    finalize_dataframe(tot, P);

    return List::create(_["windows"] = df, _["total"] = tot,
                        _["trait_names"] = ctx.trait_names,
                        _["method"] = (pr.method == 0 ? "HE" : "REML"));
}
