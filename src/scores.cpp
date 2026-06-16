// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <Rcpp.h>
#include "readBedBlock.h"
#include "geno_utils.h"
#include <vector>
#include <string>
#include <unordered_map>

using namespace Rcpp;
using Eigen::VectorXd;

// [[Rcpp::export]]
Rcpp::DataFrame compute_prs(
    const std::string& filename,       // path to bed/bim/fam (without extension)
    const Rcpp::CharacterVector& snp_input,  // SNP IDs
    const Rcpp::NumericVector&   beta_input, // corresponding beta values
    int block_size                     // SNPs per block (e.g. 256)
) {
    // ── Validate input ───────────────────────────────────────────────────────
    if (snp_input.size() != beta_input.size())
        Rcpp::stop("snp_input and beta_input must have the same length");

    // ── Read .bim and .fam ───────────────────────────────────────────────────
    List bim = read_bim_file(filename);
    List fam = read_fam_file(filename);

    Rcpp::CharacterVector snp_ids  = bim["snp"];
    Rcpp::CharacterVector geno_iid = fam["iid"];

    int n_total = geno_iid.size();
    int n_snps  = snp_ids.size();

    Rcout << "Individuals: " << n_total << "\n";
    Rcout << "SNPs in bim: " << n_snps  << "\n";

    // ── Build hash map: SNP ID -> beta ───────────────────────────────────────
    std::unordered_map<std::string, double> beta_map;
    for (int i = 0; i < snp_input.size(); ++i) {
        beta_map[Rcpp::as<std::string>(snp_input[i])] = beta_input[i];
    }
    Rcout << "SNPs provided: " << beta_map.size() << "\n";

    // ── Match beta values to bim order ───────────────────────────────────────
    std::vector<double> beta_vec(n_snps, 0.0);
    int n_matched = 0;
    for (int j = 0; j < n_snps; ++j) {
        std::string sid = Rcpp::as<std::string>(snp_ids[j]);
        auto it = beta_map.find(sid);
        if (it != beta_map.end()) {
            beta_vec[j] = it->second;
            n_matched++;
        }
    }
    Rcout << "SNPs matched between bim and input: " << n_matched << "\n";

    // ── Initialise PRS vector ────────────────────────────────────────────────
    VectorXd prs(n_total);
    prs.setZero();

    // ── Block loop ───────────────────────────────────────────────────────────
    int n_blocks = (n_snps + block_size - 1) / block_size;

    for (int b = 0; b < n_blocks; ++b) {

        int snp_start  = b * block_size;
        int snp_end    = std::min(n_snps - 1, snp_start + block_size - 1);
        int n_snps_blk = snp_end - snp_start + 1;

        // Skip block if all betas are zero
        bool any_nonzero = false;
        for (int s = 0; s < n_snps_blk; ++s) {
            if (std::abs(beta_vec[snp_start + s]) > 1e-10) {
                any_nonzero = true;
                break;
            }
        }
        if (!any_nonzero) continue;

        // Read raw genotype block
        IntegerMatrix geno_block = readBedBlock(
            filename + ".bed",
            n_total, n_snps,
            0, n_total - 1,
            snp_start, snp_end
        );

        // Loop over SNPs in block
        for (int s = 0; s < n_snps_blk; ++s) {
            int j = snp_start + s;

            double beta_j = beta_vec[j];
            if (std::abs(beta_j) < 1e-10) continue;

            // Standardise genotype column
            double gsum = 0.0, gsumsq = 0.0;
            int n_used = 0;
            for (int i = 0; i < n_total; ++i) {
                int g = geno_block(i, s);
                if (g != -1) { gsum += g; gsumsq += (double)g * g; n_used++; }
            }
            if (n_used < 3) continue;

            double gmean  = gsum / n_used;
            double gvar   = (gsumsq - n_used * gmean * gmean) / (n_used - 1);
            double inv_sd = (gvar > 1e-10) ? 1.0 / std::sqrt(gvar) : 0.0;
            if (inv_sd < 1e-10) continue;

            // Accumulate PRS
            for (int i = 0; i < n_total; ++i) {
                int g = geno_block(i, s);
                double g_std = (g != -1) ? (g - gmean) * inv_sd : 0.0;
                prs(i) += g_std * beta_j;
            }
        }

        if ((b + 1) % 100 == 0 || b == n_blocks - 1)
            Rcout << "  Block " << b + 1 << "/" << n_blocks << "\r";
    }
    Rcout << "\nDone.\n";

    // ── Return data frame with IID and PRS ───────────────────────────────────
    return Rcpp::DataFrame::create(
        Rcpp::Named("IID") = geno_iid,
        Rcpp::Named("PRS") = Rcpp::wrap(prs)
    );
}