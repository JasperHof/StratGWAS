// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppParallel)]]
#include <RcppEigen.h>
#include <Rcpp.h>
#include <RcppParallel.h>
#include "readBedBlock.h"
#include "geno_utils.h"
#include <vector>
#include <algorithm>
#include <cmath>

using namespace Rcpp;
using namespace RcppParallel;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// Soft-threshold for lasso component
inline double soft_threshold(double x, double lambda) {
    if (x >  lambda) return x - lambda;
    if (x < -lambda) return x + lambda;
    return 0.0;
}

// Single-SNP VB update
inline double update_snp_vb(
    double  XtR_j,
    double  XtX_j,
    double& mu_j,
    double& sigma2_j,
    double  h2_j,
    double  p_en,
    double  F_en,
    double  sigma2_e
) {
    double mu_old = mu_j;

    double ridge_var  = (F_en > 1e-10 && (1.0 - p_en) > 1e-10)
                        ? F_en * h2_j / (1.0 - p_en) : 1e-10;
    double lasso_rate = (p_en > 1e-10 && (1.0 - F_en) > 1e-10 && h2_j > 1e-10)
                        ? std::sqrt(2.0 * p_en / ((1.0 - F_en) * h2_j)) : 1e8;

    if (p_en < 1e-8) {
        sigma2_j = 1.0 / (XtX_j / sigma2_e + 1.0 / ridge_var);
        mu_j     = sigma2_j * XtR_j / sigma2_e;

    } else if (p_en > 1.0 - 1e-8) {
        double naive = XtR_j / XtX_j;
        mu_j         = soft_threshold(naive, sigma2_e * lasso_rate / XtX_j);
        sigma2_j     = sigma2_e / XtX_j;

    } else {
        double sigma2_ridge = 1.0 / (XtX_j / sigma2_e + 1.0 / ridge_var);
        double mu_ridge     = sigma2_ridge * XtR_j / sigma2_e;

        double sigma2_lasso = sigma2_e / XtX_j;
        double mu_lasso     = soft_threshold(XtR_j / XtX_j,
                                             sigma2_e * lasso_rate / XtX_j);

        double log_w_ridge = std::log(1.0 - p_en)
                           + 0.5 * mu_ridge * mu_ridge / sigma2_ridge
                           - 0.5 * std::log(ridge_var)
                           + 0.5 * std::log(sigma2_ridge);

        double log_w_lasso = std::log(p_en)
                           + std::log(lasso_rate / 2.0)
                           + 0.5 * mu_lasso * mu_lasso / sigma2_lasso;

        double log_max = std::max(log_w_ridge, log_w_lasso);
        double w_ridge = std::exp(log_w_ridge - log_max);
        double w_lasso = std::exp(log_w_lasso - log_max);
        double Z       = w_ridge + w_lasso;
        w_ridge /= Z;
        w_lasso /= Z;

        mu_j     = w_ridge * mu_ridge + w_lasso * mu_lasso;
        sigma2_j = w_ridge * (sigma2_ridge + mu_ridge * mu_ridge)
                 + w_lasso * (sigma2_lasso + mu_lasso * mu_lasso)
                 - mu_j * mu_j;
        sigma2_j = std::max(sigma2_j, 1e-10);
    }

    return std::abs(mu_j - mu_old);
}

// Worker: parallel VB updates over phenotypes for one genotype chunk.
//
// RMatrix / RVector only wrap Rcpp types (NumericMatrix etc.), NOT Eigen types.
// So we store Eigen data as raw const pointers and access via pointer arithmetic.
// Mutable per-phenotype state (mu_blk, sigma2_blk, pgs_mat) is stored as
// RMatrix<double> wrapping NumericMatrix -- exactly like PhenotypeWorker in
// your GWAS code. Each thread owns a disjoint column range so no locking needed.

struct VBChunkWorker : public Worker {

    // Read-only: raw pointers into Eigen storage (column-major)
    const double* G_ptr;        // n_inds x n_snps_blk
    const double* XtX_ptr;      // n_snps_blk
    const double* Y_ptr;        // n_inds x n_pheno
    const int n_inds;
    const int n_snps_blk;
    const int n_pheno;

    // Read-only scalar/vector parameters
    const std::vector<double>& h2_j_blk;
    const std::vector<int>&    chr_blk;
    const double p_en;
    const double F_en;
    const double sigma2_e;
    const bool   loco;
    const int    loco_chr;

    // Mutable per-phenotype state: RMatrix wraps NumericMatrix (Rcpp type - OK)
    RMatrix<double> mu_blk;
    RMatrix<double> sigma2_blk;
    RMatrix<double> pgs_mat;

    // Output: chunk LL per phenotype
    std::vector<double>& chunk_ll;

    VBChunkWorker(
        const MatrixXd&            G,
        const VectorXd&            XtX_vec,
        const MatrixXd&            Y_mat,
        int n_inds_, int n_snps_blk_, int n_pheno_,
        const std::vector<double>& h2_j_blk_,
        const std::vector<int>&    chr_blk_,
        double p_en_, double F_en_, double sigma2_e_,
        bool loco_, int loco_chr_,
        NumericMatrix&             mu_blk_nm,
        NumericMatrix&             sigma2_blk_nm,
        NumericMatrix&             pgs_nm,
        std::vector<double>&       chunk_ll_
    )
      : G_ptr(G.data()), XtX_ptr(XtX_vec.data()), Y_ptr(Y_mat.data()),
        n_inds(n_inds_), n_snps_blk(n_snps_blk_), n_pheno(n_pheno_),
        h2_j_blk(h2_j_blk_), chr_blk(chr_blk_),
        p_en(p_en_), F_en(F_en_), sigma2_e(sigma2_e_),
        loco(loco_), loco_chr(loco_chr_),
        mu_blk(mu_blk_nm), sigma2_blk(sigma2_blk_nm), pgs_mat(pgs_nm),
        chunk_ll(chunk_ll_) {}

    void operator()(std::size_t begin, std::size_t end) {

        for (std::size_t p = begin; p < end; ++p) {

            double ll = 0.0;

            // Column pointers for this phenotype
            // Eigen MatrixXd is column-major: column p starts at offset p*n_rows
            const double* Y_p   = Y_ptr   + p * n_inds;
            double*       pgs_p = &pgs_mat(0, p);      // RMatrix column p

            for (int s = 0; s < n_snps_blk; ++s) {

                if (loco && chr_blk[s] == loco_chr) continue;
                if (h2_j_blk[s] <= 0.0)             continue;

                double XtX_j = XtX_ptr[s];
                if (XtX_j < 1e-6) continue;

                // G column s: offset s*n_inds in column-major storage
                const double* G_s = G_ptr + s * n_inds;

                double mu_sp = mu_blk(s, p);

                // Compute X_s^T * residual, where residual = Y_p - pgs_p + X_s*mu_sp
                double XtR_j = 0.0;
                for (int i = 0; i < n_inds; ++i) {
                    double r = Y_p[i] - pgs_p[i] + G_s[i] * mu_sp;
                    XtR_j += G_s[i] * r;
                }

                double mu_old    = mu_sp;
                double sigma2_sp = sigma2_blk(s, p);

                double delta = update_snp_vb(
                    XtR_j, XtX_j,
                    mu_blk(s, p), sigma2_sp,
                    h2_j_blk[s], p_en, F_en, sigma2_e
                );
                sigma2_blk(s, p) = sigma2_sp;

                // Update PGS for phenotype p
                double d_mu = mu_blk(s, p) - mu_old;
                for (int i = 0; i < n_inds; ++i)
                    pgs_p[i] += G_s[i] * d_mu;

                ll += delta;
            }

            chunk_ll[p] = ll;
        }
    }
};

// [[Rcpp::export]]
Rcpp::List vb_elastic_net_prs_multi(
    const std::string& filename,
    const SEXP         pheno_mat,
    const Rcpp::NumericVector& maf_all,
    const Rcpp::IntegerVector& chr_all,
    double h2,
    double alpha_param,
    double p_en,
    double F_en,
    int    chunk_size,
    double tol,
    int    max_scans,
    bool   loco,
    int    loco_chr
) {
    List bim = read_bim_file(filename);
    List fam = read_fam_file(filename);

    Rcpp::CharacterVector snp     = bim["snp"];
    CharacterVector geno_iid      = fam["iid"];

    int n_total_inds = geno_iid.size();
    int n_snps       = snp.size();

    Rcout << "Individuals (from .fam): " << n_total_inds << "\n";
    Rcout << "SNPs (from .bim): "        << n_snps       << "\n";

    if (!Rf_isMatrix(pheno_mat) || Rf_isNull(rownames(pheno_mat)))
        stop("Phenotype must be a numeric matrix with IIDs as rownames");

    NumericMatrix pheno    = as<NumericMatrix>(pheno_mat);
    CharacterVector pheno_ids = rownames(pheno);
    int n_pheno = pheno.ncol();

    Rcout << "Phenotypes: " << n_pheno << "\n";

    // Match individuals
    IntegerVector match_idx = match(geno_iid, pheno_ids);
    std::vector<int> geno_keep, pheno_keep;
    for (int i = 0; i < match_idx.size(); ++i) {
        if (match_idx[i] != NA_INTEGER) {
            geno_keep.push_back(i);
            pheno_keep.push_back(match_idx[i] - 1);
        }
    }
    if (geno_keep.empty())
        stop("No overlapping individuals between genotype and phenotype");

    int n_inds = (int)geno_keep.size();
    Rcout << "Using " << n_inds << " individuals\n";

    // Build standardised phenotype matrix Y (n_inds x n_pheno, column-major Eigen)
    MatrixXd Y_mat(n_inds, n_pheno);
    for (int p = 0; p < n_pheno; ++p) {
        double s = 0.0;
        for (int i = 0; i < n_inds; ++i) {
            Y_mat(i, p) = pheno(pheno_keep[i], p);
            s += Y_mat(i, p);
        }
        double mean = s / n_inds;
        Y_mat.col(p).array() -= mean;
        double sd = std::sqrt(Y_mat.col(p).squaredNorm() / (n_inds - 1));
        if (sd > 1e-10) Y_mat.col(p) /= sd;
    }

    // Per-SNP heritability weights
    std::vector<double> h2_j(n_snps, 0.0);
    double W = 0.0;
    for (int j = 0; j < n_snps; ++j) {
        if (loco && chr_all[j] == loco_chr) continue;
        double f = maf_all[j];
        if (f <= 0.0 || f >= 1.0) continue;
        double w = std::pow(f * (1.0 - f), 1.0 + alpha_param);
        h2_j[j] = w;
        W += w;
    }
    if (W > 0.0)
        for (int j = 0; j < n_snps; ++j) h2_j[j] *= h2 / W;

    double sigma2_e = std::max(1.0 - h2, 1e-6);

    // Global VB state as NumericMatrix (Rcpp type so RMatrix can wrap them)
    NumericMatrix mu_global(n_snps, n_pheno);
    NumericMatrix sigma2_global(n_snps, n_pheno);
    NumericMatrix pgs_global(n_inds, n_pheno);
    std::fill(mu_global.begin(),     mu_global.end(),     0.0);
    std::fill(sigma2_global.begin(), sigma2_global.end(), 1e-6);
    std::fill(pgs_global.begin(),    pgs_global.end(),    0.0);

    // Chunk bookkeeping
    int n_chunks = (n_snps + chunk_size - 1) / chunk_size;
    std::vector<bool>   chunk_active(n_chunks * n_pheno, true);
    std::vector<double> chunk_ll_prev(n_chunks * n_pheno, 1e30);

    long long total_updates = 0;

    for (int scan = 0; scan < max_scans; ++scan) {

        bool any_active = false;
        for (int cp = 0; cp < n_chunks * n_pheno; ++cp)
            if (chunk_active[cp]) { any_active = true; break; }
        if (!any_active) break;

        int n_active_chunks = 0;

        for (int c = 0; c < n_chunks; ++c) {

            bool chunk_needed = false;
            for (int p = 0; p < n_pheno; ++p)
                if (chunk_active[c * n_pheno + p]) { chunk_needed = true; break; }
            if (!chunk_needed) continue;
            ++n_active_chunks;

            int snp_start  = c * chunk_size;
            int snp_end    = std::min(n_snps - 1, snp_start + chunk_size - 1);
            int n_snps_blk = snp_end - snp_start + 1;

            // Read genotype block once -- shared across all phenotypes
            IntegerMatrix geno_block_full = readBedBlock(
                filename + ".bed",
                n_total_inds, n_snps,
                0, n_total_inds - 1,
                snp_start, snp_end
            );

            IntegerMatrix geno_sub(n_inds, n_snps_blk);
            for (int i = 0; i < n_inds; ++i)
                for (int s = 0; s < n_snps_blk; ++s)
                    geno_sub(i, s) = geno_block_full(geno_keep[i], s);

            // Standardise genotypes -- shared
            MatrixXd G(n_inds, n_snps_blk);
            for (int s = 0; s < n_snps_blk; ++s) {
                double gsum = 0.0, gsumsq = 0.0;
                int n_used = 0;
                for (int i = 0; i < n_inds; ++i) {
                    int g = geno_sub(i, s);
                    if (g != -1) { gsum += g; gsumsq += (double)g * g; n_used++; }
                }
                if (n_used < 3) { G.col(s).setZero(); continue; }
                double gmean  = gsum / n_used;
                double gvar   = (gsumsq - n_used * gmean * gmean) / (n_used - 1);
                double inv_sd = (gvar > 1e-10) ? 1.0 / std::sqrt(gvar) : 0.0;
                for (int i = 0; i < n_inds; ++i) {
                    int g = geno_sub(i, s);
                    G(i, s) = (g != -1) ? (g - gmean) * inv_sd : 0.0;
                }
            }

            // Precompute XtX -- shared
            VectorXd XtX_vec(n_snps_blk);
            for (int s = 0; s < n_snps_blk; ++s)
                XtX_vec(s) = G.col(s).squaredNorm();

            // Block-level heritability and chromosome
            std::vector<double> h2_j_blk(n_snps_blk);
            std::vector<int>    chr_blk(n_snps_blk);
            for (int s = 0; s < n_snps_blk; ++s) {
                h2_j_blk[s] = h2_j[snp_start + s];
                chr_blk[s]  = chr_all[snp_start + s];
            }

            // Slice global mu / sigma2 into a block-sized NumericMatrix
            // so the Worker indexes with (s, p) not (snp_start+s, p)
            NumericMatrix mu_blk(n_snps_blk, n_pheno);
            NumericMatrix sigma2_blk(n_snps_blk, n_pheno);
            for (int s = 0; s < n_snps_blk; ++s)
                for (int p = 0; p < n_pheno; ++p) {
                    mu_blk(s, p)     = mu_global(snp_start + s, p);
                    sigma2_blk(s, p) = sigma2_global(snp_start + s, p);
                }

            std::vector<double> chunk_ll_out(n_pheno, 0.0);

            // Parallel VB updates over phenotypes
            VBChunkWorker worker(
                G, XtX_vec, Y_mat,
                n_inds, n_snps_blk, n_pheno,
                h2_j_blk, chr_blk,
                p_en, F_en, sigma2_e,
                loco, loco_chr,
                mu_blk, sigma2_blk, pgs_global,
                chunk_ll_out
            );
            parallelFor(0, n_pheno, worker);

            // Write block mu / sigma2 back to global storage
            for (int s = 0; s < n_snps_blk; ++s)
                for (int p = 0; p < n_pheno; ++p) {
                    mu_global(snp_start + s, p)     = mu_blk(s, p);
                    sigma2_global(snp_start + s, p) = sigma2_blk(s, p);
                }

            // Update convergence flags
            for (int p = 0; p < n_pheno; ++p) {
                int cp = c * n_pheno + p;
                double ll_change = std::abs(chunk_ll_out[p] - chunk_ll_prev[cp]);
                chunk_ll_prev[cp] = chunk_ll_out[p];
                if (ll_change < (double)n_inds * tol)
                    chunk_active[cp] = false;
                total_updates += n_snps_blk;
            }
        }

        int still_active = 0;
        for (int cp = 0; cp < n_chunks * n_pheno; ++cp)
            if (chunk_active[cp]) ++still_active;

        Rcout << "Scan " << scan + 1
              << ": " << n_active_chunks << "/" << n_chunks
              << " active chunks | "
              << still_active << "/" << (n_chunks * n_pheno)
              << " (chunk x pheno) pairs remaining | "
              << total_updates << " total updates\n";
    }

    return Rcpp::List::create(
        Rcpp::Named("pgs")       = pgs_global,
        Rcpp::Named("beta")      = mu_global,
        Rcpp::Named("beta_var")  = sigma2_global,
        Rcpp::Named("snp")       = snp,
        Rcpp::Named("n_updates") = (double)total_updates
    );
}