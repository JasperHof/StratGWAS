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

using namespace Rcpp;

using Eigen::Matrix2d;
using Eigen::Vector2d;

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

    // Read in phenotype
    Rcpp::NumericMatrix pheno;
    CharacterVector pheno_ids;

    if (Rf_isMatrix(pheno_mat) && !Rf_isNull(rownames(pheno_mat))) { // Case 1: rownames as IDs
        pheno = as<NumericMatrix>(pheno_mat);
        pheno_ids = rownames(pheno_mat);
    } else {
        stop("Phenotype must be a numerical matrix with IDs as rownames");
    }

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
    
    // Reorder phenotype to match genotype
    NumericMatrix pheno_sub(geno_keep.size(), n_pheno);
    for (size_t i = 0; i < geno_keep.size(); ++i) {
        for (int j = 0; j < n_pheno; ++j) {
            pheno_sub(i, j) = pheno(pheno_keep[i], j);
        }
    }

    pheno = pheno_sub;

    int n_inds = pheno.nrow(); // note that this is different from n_total_inds (genotype dim)
    Rcout << "Will be using " << n_inds << " overlapping individuals\n";

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

        // Standardize genotypes
        for (int s = 0; s < G.cols(); ++s) {
            Eigen::VectorXd v = G.col(s);
            double mean = v.mean();
            double sd = std::sqrt((v.array() - mean).square().sum() / (v.size() - 1));
            if (sd > 0)
                G.col(s) = (v.array() - mean) / sd;
            else
                G.col(s).setZero();
        }

        // Perform the linear regression for each phenotype
        for (int p = 0; p < n_pheno; ++p) {
                
            Eigen::VectorXd y = Eigen::Map<Eigen::VectorXd>(
                pheno.column(p).begin(), n_inds);

            // Pre-allocate SNP
            std::vector<double> xbuf(n_inds);
            std::vector<double> ybuf(n_inds);

            // now loop over the SNPs
            for (int s = 0; s < n_snps_block; ++s) {

                int n_used = 0;

                for (int i = 0; i < n_inds; ++i) {
                    if (geno_sub(i, s) != -1 && !NumericVector::is_na(y(i))) {
                        xbuf[n_used] = G(i, s);
                        ybuf[n_used] = y(i);
                        n_used++;
                    }
                }

                double beta = NA_REAL, se = NA_REAL, chisq = NA_REAL, pval = NA_REAL;

                if (n_used > 2) {
                    Eigen::Map<Eigen::VectorXd> xs(xbuf.data(), n_used);
                    Eigen::Map<Eigen::VectorXd> ys(ybuf.data(), n_used);

                    double xTx = xs.dot(xs);
                    if (xTx > 0) {
                        beta = xs.dot(ys) / xTx;
                        double rss = (ys - xs * beta).squaredNorm();
                        double sigma2 = rss / (n_used - 2);
                        se = std::sqrt(sigma2 / xTx);
                        chisq = (beta / se) * (beta / se);
                        pval = R::pchisq(chisq, 1, false, false);
                    }
                }

                out_files[p]
                    << chr[block_start + s] << "\t"
                    << snp[block_start + s] << "\t"
                    << pos[block_start + s] << "\t"
                    << a1[block_start + s]  << "\t"
                    << a2[block_start + s]  << "\t"
                    << beta << "\t"
                    << se << "\t"
                    << chisq << "\t"
                    << pval << "\t"
                    << n_used << "\t"
                    << maf[s] << "\t"
                    << miss[s] << "\n";
            }
        }
    }

    return 0;
}