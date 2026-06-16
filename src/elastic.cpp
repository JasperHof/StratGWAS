// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <Rcpp.h>
#include "readBedBlock.h"
#include "geno_utils.h"
#include <vector>
#include <random>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>

using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// ─────────────────────────────────────────────────────────────────────────────
// Soft-threshold for lasso component
// ─────────────────────────────────────────────────────────────────────────────
inline double soft_threshold(double x, double lambda) {
    if (x >  lambda) return x - lambda;
    if (x < -lambda) return x + lambda;
    return 0.0;
}

// ─────────────────────────────────────────────────────────────────────────────
// Single-SNP variational Bayes update under elastic net prior
// Updates mu_j and sigma2_j in place, returns |delta mu| for convergence
// ─────────────────────────────────────────────────────────────────────────────
inline double update_snp_vb(
    double  XtR_j,     // X_j^T * residual (with current SNP added back)
    double  XtX_j,     // ||X_j||^2
    double& mu_j,      // posterior mean (in/out)
    double& sigma2_j,  // posterior variance (in/out)
    double  h2_j,      // per-SNP heritability
    double  p_en,      // lasso proportion
    double  F_en,      // ridge proportion
    double  sigma2_e   // residual variance
) {
    double mu_old = mu_j;

    // Prior variances from Equation (4) of LDAK-KVIK supplemental
    double ridge_var  = (F_en > 1e-10 && (1.0 - p_en) > 1e-10)
                        ? F_en * h2_j / (1.0 - p_en) : 1e-10;
    double lasso_rate = (p_en > 1e-10 && (1.0 - F_en) > 1e-10 && h2_j > 1e-10)
                        ? std::sqrt(2.0 * p_en / ((1.0 - F_en) * h2_j)) : 1e8;

    if (p_en < 1e-8) {
        // ── Pure ridge ──────────────────────────────────────────────────────
        sigma2_j = 1.0 / (XtX_j / sigma2_e + 1.0 / ridge_var);
        mu_j     = sigma2_j * XtR_j / sigma2_e;

    } else if (p_en > 1.0 - 1e-8) {
        // ── Pure lasso ──────────────────────────────────────────────────────
        double naive = XtR_j / XtX_j;
        mu_j         = soft_threshold(naive, sigma2_e * lasso_rate / XtX_j);
        sigma2_j     = sigma2_e / XtX_j;

    } else {
        // ── Elastic net mixture ─────────────────────────────────────────────

        // Ridge component posterior
        double sigma2_ridge = 1.0 / (XtX_j / sigma2_e + 1.0 / ridge_var);
        double mu_ridge     = sigma2_ridge * XtR_j / sigma2_e;

        // Lasso component posterior (Gaussian approximation to DE)
        double sigma2_lasso = sigma2_e / XtX_j;
        double mu_lasso     = soft_threshold(XtR_j / XtX_j,
                                             sigma2_e * lasso_rate / XtX_j);

        // Log unnormalised mixture weights
        double log_w_ridge = std::log(1.0 - p_en)
                           + 0.5 * mu_ridge * mu_ridge / sigma2_ridge
                           - 0.5 * std::log(ridge_var)
                           + 0.5 * std::log(sigma2_ridge);

        double log_w_lasso = std::log(p_en)
                           + std::log(lasso_rate / 2.0)
                           + 0.5 * mu_lasso * mu_lasso / sigma2_lasso;

        // Numerically stable softmax
        double log_max = std::max(log_w_ridge, log_w_lasso);
        double w_ridge = std::exp(log_w_ridge - log_max);
        double w_lasso = std::exp(log_w_lasso - log_max);
        double Z       = w_ridge + w_lasso;
        w_ridge /= Z;
        w_lasso /= Z;

        // Combined posterior moments (law of total expectation/variance)
        mu_j     = w_ridge * mu_ridge + w_lasso * mu_lasso;
        sigma2_j = w_ridge * (sigma2_ridge + mu_ridge * mu_ridge)
                 + w_lasso * (sigma2_lasso + mu_lasso * mu_lasso)
                 - mu_j * mu_j;
        sigma2_j = std::max(sigma2_j, 1e-10);
    }

    return std::abs(mu_j - mu_old);
}

// [[Rcpp::export]]
Rcpp::List vb_elastic_net_prs(
    const std::string& filename,
    const SEXP         pheno_mat,    // numeric matrix, rownames = IIDs
    const Rcpp::NumericVector& maf_all,  // per-SNP MAF, length n_snps
    const Rcpp::IntegerVector& chr_all,  // per-SNP chromosome, length n_snps
    double h2,           // estimated SNP heritability
    double alpha_param,  // MAF-weighting exponent (e.g. -0.25)
    double p_en,         // elastic net lasso proportion
    double F_en,         // elastic net ridge proportion
    int    chunk_size,   // SNPs per chunk (e.g. 256)
    double tol,          // convergence: chunk LL change < n * tol
    int    max_scans,    // maximum genome scans
    bool   loco,         // LOCO mode?
    int    loco_chr      // chromosome to exclude if loco = true
) {
    // ── Read .bim and .fam ───────────────────────────────────────────────────
    List bim = read_bim_file(filename);
    List fam = read_fam_file(filename);

    Rcpp::CharacterVector snp     = bim["snp"];
    Rcpp::IntegerVector   chr_bim = bim["chr"];

    CharacterVector geno_iid    = fam["iid"];
    int n_total_inds = geno_iid.size();
    int n_snps       = snp.size();

    Rcout << "Number of individuals (from .fam): " << n_total_inds << "\n";
    Rcout << "Number of SNPs (from .bim): "        << n_snps       << "\n";

    // ── Read and validate phenotype ──────────────────────────────────────────
    if (!Rf_isMatrix(pheno_mat) || Rf_isNull(rownames(pheno_mat)))
        stop("Phenotype must be a numeric matrix with IIDs as rownames");

    NumericMatrix pheno    = as<NumericMatrix>(pheno_mat);
    CharacterVector pheno_ids = rownames(pheno);

    // ── Match individuals ────────────────────────────────────────────────────
    IntegerVector match_idx = match(geno_iid, pheno_ids);
    std::vector<int> geno_keep, pheno_keep;
    for (int i = 0; i < match_idx.size(); ++i) {
        if (match_idx[i] != NA_INTEGER) {
            geno_keep.push_back(i);
            pheno_keep.push_back(match_idx[i] - 1);
        }
    }
    if (geno_keep.empty()) stop("No overlapping individuals between genotype and phenotype");

    int n_inds = (int)geno_keep.size();
    Rcout << "Using " << n_inds << " individuals\n";

    // ── Build phenotype vector (first column only), centre and scale ─────────
    VectorXd Y(n_inds);
    for (int i = 0; i < n_inds; ++i) Y(i) = pheno(pheno_keep[i], 0);
    Y.array() -= Y.mean();
    double y_sd = std::sqrt(Y.squaredNorm() / (n_inds - 1));
    if (y_sd > 1e-10) Y /= y_sd;

    // ── Per-SNP heritability weights: h2_j = w_j * h2 / W ───────────────────
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

    // ── Initialise variational parameters and PGS ────────────────────────────
    VectorXd mu(n_snps);      mu.setZero();
    VectorXd sigma2(n_snps);  sigma2.setConstant(1e-6);
    VectorXd pgs(n_inds);     pgs.setZero();

    // ── Chunk bookkeeping ─────────────────────────────────────────────────────
    int n_chunks = (n_snps + chunk_size - 1) / chunk_size;
    std::vector<bool>   chunk_active(n_chunks, true);
    std::vector<double> chunk_ll_prev(n_chunks, 1e30);
    int total_updates = 0;

    // ── Outer scan loop ───────────────────────────────────────────────────────
    for (int scan = 0; scan < max_scans; ++scan) {

        int  n_active   = 0;
        bool any_active = false;

        for (int c = 0; c < n_chunks; ++c) {
            if (!chunk_active[c]) continue;
            any_active = true;
            n_active++;

            int snp_start  = c * chunk_size;
            int snp_end    = std::min(n_snps - 1, snp_start + chunk_size - 1);
            int n_snps_blk = snp_end - snp_start + 1;

            // ── Read raw genotype block ──────────────────────────────────────
            IntegerMatrix geno_block = readBedBlock(
                filename + ".bed",
                n_total_inds, n_snps,
                0, n_total_inds - 1,
                snp_start, snp_end
            );

            // ── Subset to matched individuals ────────────────────────────────
            IntegerMatrix geno_sub(n_inds, n_snps_blk);
            for (int i = 0; i < n_inds; ++i)
                for (int s = 0; s < n_snps_blk; ++s)
                    geno_sub(i, s) = geno_block(geno_keep[i], s);

            // ── Standardise each SNP column ──────────────────────────────────
            MatrixXd G(n_inds, n_snps_blk);
            Rcpp::NumericVector maf_blk  = computeMAF(geno_sub);
            Rcpp::NumericVector miss_blk = computeMissingness(geno_sub);

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

            // ── VB updates for each SNP in chunk ─────────────────────────────
            double chunk_ll = 0.0;

            for (int s = 0; s < n_snps_blk; ++s) {
                int j = snp_start + s;

                if (loco && chr_all[j] == loco_chr) continue;
                if (h2_j[j] <= 0.0) continue;

                double XtX_j = G.col(s).squaredNorm();
                if (XtX_j < 1e-6) continue;

                // Add current SNP back into residual before update
                VectorXd R = Y - pgs + G.col(s) * mu(j);

                double XtR_j = G.col(s).dot(R);
                double mu_old = mu(j);

                // VB update
                double delta = update_snp_vb(
                    XtR_j, XtX_j,
                    mu(j), sigma2(j),
                    h2_j[j], p_en, F_en, sigma2_e
                );

                // Update PGS
                pgs.noalias() += G.col(s) * (mu(j) - mu_old);

                chunk_ll += delta;
                total_updates++;
            }

            // ── Chunk convergence check ──────────────────────────────────────
            double ll_change = std::abs(chunk_ll - chunk_ll_prev[c]);
            chunk_ll_prev[c] = chunk_ll;
            if (ll_change < (double)n_inds * tol)
                chunk_active[c] = false;
        }

        Rcout << "Scan " << scan + 1 << ": "
              << n_active << "/" << n_chunks << " active chunks, "
              << total_updates << " total updates\n";

        if (!any_active) break;
    }

    // ── Return ────────────────────────────────────────────────────────────────
    return Rcpp::List::create(
        Rcpp::Named("pgs")       = Rcpp::wrap(pgs),
        Rcpp::Named("beta")      = Rcpp::wrap(mu),
        Rcpp::Named("beta_var")  = Rcpp::wrap(sigma2),
        Rcpp::Named("snp")       = snp,
        Rcpp::Named("n_updates") = total_updates
    );
}
