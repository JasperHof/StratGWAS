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
using Eigen::Matrix2d;
using Eigen::Vector2d;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// [[Rcpp::export]]
Rcpp::NumericMatrix linear_gwas(const std::string& filename, const SEXP pheno_mat, int block_size, const std::string& out_file) {
    
    // Read in .bim and .fam file
    List bim = read_bim_file(filename);
    List fam = read_fam_file(filename);

    Rcpp::CharacterVector snp = bim["snp"];
    Rcpp::IntegerVector chr = bim["chr"];
    Rcpp::IntegerVector pos = bim["pos"];
    Rcpp::CharacterVector a1 = bim["a1"];
    Rcpp::CharacterVector a2 = bim["a2"];

    CharacterVector geno_iid = fam["iid"];
    int n_total_inds = geno_iid.size();
    int n_snps = snp.size();

    Rcpp::Rcout << "Number of individuals (from .fam): " << n_total_inds << std::endl;
    Rcpp::Rcout << "Number of SNPs (from .bim): " << n_snps << std::endl;   

    // Read phenotype
    if (!Rf_isMatrix(pheno_mat) || Rf_isNull(rownames(pheno_mat)))
        stop("Phenotype must be a numerical matrix with IDs as rownames");

    NumericMatrix pheno = as<NumericMatrix>(pheno_mat);
    CharacterVector pheno_ids = rownames(pheno);
    int n_pheno = pheno.cols();

    // Match phenotype ids to genotype ids
    IntegerVector match_idx = match(geno_iid, pheno_ids);

    std::vector<int> geno_keep;
    std::vector<int> pheno_keep;

    for (int i = 0; i < match_idx.size(); ++i) {
        if (match_idx[i] != NA_INTEGER) {
            geno_keep.push_back(i);
            pheno_keep.push_back(match_idx[i] - 1); // R → C++ index
        }
    }

    if (geno_keep.empty())
        stop("No overlapping individuals between genotype and phenotype");

    int n_inds = geno_keep.size();
    Rcout << "Using " << n_inds << " individuals\n";

    // Reorder phenotype to match genotype
    NumericMatrix pheno_sub(n_inds, n_pheno);
    for (size_t i = 0; i < n_inds; ++i) {
        for (int j = 0; j < n_pheno; ++j) {
            pheno_sub(i, j) = pheno(pheno_keep[i], j);
        }
    }

    pheno = pheno_sub;

    // Open output files for each phenotype and write header
    std::vector<std::ofstream> out_files(n_pheno);
    for(int i = 0; i < n_pheno; ++i) {
        std::string fname = out_file + ".pheno" + std::to_string(i+1);
        out_files[i].open(fname, std::ios::out);
        out_files[i] << "Chromosome\tPredictor\tBasepair\tA1\tA2\tBeta\tSE\tChisq\tP\tN\tMAF\tMiss\n";
    }

    int nr_blocks = n_snps / block_size;

    for(int b = 0; b < nr_blocks; ++b){

        Rcout << "Performing linear regression for block " << b + 1 << "/" << nr_blocks << "\n";

        int block_start = b * block_size;
        int block_end = std::min(n_snps - 1, (b + 1) * block_size - 1); // make sure you don't exceed n_snps
        int n_snps_block = block_end - block_start + 1;

        // Read in the genotype and subset
        IntegerMatrix geno_block = readBedBlock(
            filename + ".bed",
            n_total_inds, n_snps,
            0, n_total_inds - 1,
            block_start, block_end
        );

        IntegerMatrix geno_sub(n_inds, n_snps_block);
        for (int i = 0; i < n_inds; ++i)
            for (int s = 0; s < n_snps_block; ++s)
                geno_sub(i, s) = geno_block(geno_keep[i], s);
        
        Eigen::MatrixXd G = as<Eigen::MatrixXd>(geno_sub);

        // Compute MAF and missingness
        Rcpp::NumericVector maf = computeMAF(geno_sub);
        Rcpp::NumericVector miss = computeMissingness(geno_sub);

        // Perform the linear regression for each phenotype
        for (int p = 0; p < n_pheno; ++p) {
                
            // Phenotype vector + mask
            VectorXd y(n_inds);
            VectorXd y_mask(n_inds);

            int n_y = 0;
            double y_sum = 0.0;

            for (int i = 0; i < n_inds; ++i) {
                if (NumericVector::is_na(pheno(i, p))) {
                    y(i) = 0.0;
                    y_mask(i) = 0.0;
                } else {
                    y(i) = pheno(i, p);
                    y_mask(i) = 1.0;
                    y_sum += y(i);
                    n_y++;
                }
            }

            if (n_y < 3) continue;

            double y_mean = y_sum / n_y;
            double yTy = 0.0;

            for (int i = 0; i < n_inds; ++i) {
                if (y_mask(i)) {
                    y(i) -= y_mean;
                    yTy += y(i) * y(i);
                }
            }

            // -----------------------------
            // Standardize genotype block
            // -----------------------------

            MatrixXd G(n_inds, n_snps_block);
            std::vector<int> n_used_snp(n_snps_block);

            for (int s = 0; s < n_snps_block; ++s) {

                double sum = 0.0, sumsq = 0.0;
                int n_used = 0;

                for (int i = 0; i < n_inds; ++i) {
                    int g = geno_sub(i, s);
                    if (g != -1 && y_mask(i)) {
                        sum += g;
                        sumsq += g * g;
                        n_used++;
                    }
                }

                n_used_snp[s] = n_used;
                if (n_used < 3) {
                    G.col(s).setZero();
                    continue;
                }

                double mean = sum / n_used;
                double var = (sumsq - n_used * mean * mean) / (n_used - 1);
                double inv_sd = (var > 0) ? 1.0 / std::sqrt(var) : 0.0;

                for (int i = 0; i < n_inds; ++i) {
                    int g = geno_sub(i, s);
                    if (g != -1 && y_mask(i))
                        G(i, s) = (g - mean) * inv_sd;
                    else
                        G(i, s) = 0.0;
                }
            }

            // Vectorized regression
            VectorXd XtY = G.transpose() * y;
            std::ostringstream buffer;

            // buffer.reserve(1 << 20);

            for (int s = 0; s < n_snps_block; ++s) {

                int n_used = n_used_snp[s];
                if (n_used < 3) continue;

                double XtX = G.col(s).squaredNorm();
                if (XtX <= 0) continue;

                double beta = XtY(s) / XtX;
                double rss = yTy - beta * XtY(s);
                if (rss <= 0) continue;

                double se = std::sqrt((rss / (n_used - 2)) / XtX);
                double chisq = (beta / se) * (beta / se);
                double pval = R::pchisq(chisq, 1, false, false);

                int idx = block_start + s;

                buffer
                    << chr[idx] << '\t'
                    << snp[idx] << '\t'
                    << pos[idx] << '\t'
                    << a1[idx] << '\t'
                    << a2[idx] << '\t'
                    << beta << '\t'
                    << se << '\t'
                    << chisq << '\t'
                    << pval << '\t'
                    << n_used << '\t'
                    << maf[s] << '\t'
                    << miss[s] << '\n';
            }

            out_files[p] << buffer.str();
        }
    }

    return 0;
}