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
Rcpp::NumericMatrix linear_gwas(const std::string& filename, const std::string& pheno_file, const std::string& out_file, int n_inds, int n_snps, int block_size) {
    
    std::string bed_file = filename + ".bed";
    std::string bim_file = filename + ".bim";
    std::string fam_file = filename + ".fam";

    // Read phenotype file using R's read.table and make it numerical matrix
    DataFrame pheno_df = Rcpp::as<DataFrame>(Rcpp::Function("read.table")(pheno_file, Named("header", true)));
    
    int n_pheno = pheno_df.size() - 2;       // number of columns
    int nr_blocks = n_snps / block_size; // note that we actually want to have floor, which is automatically done

    NumericMatrix pheno(n_inds, n_pheno);
    for (int j = 0; j < n_pheno; ++j) {
        Rcpp::NumericVector col = pheno_df[j + 2];
        pheno(_, j) = col;
    }

    // Read BIM file into Rcpp vectors
    Rcpp::List bim = read_bim_file(filename);
    Rcpp::CharacterVector snp = bim["snp"];
    Rcpp::IntegerVector chr = bim["chr"];
    Rcpp::IntegerVector pos = bim["pos"];
    Rcpp::CharacterVector a1 = bim["a1"];
    Rcpp::CharacterVector a2 = bim["a2"];

    // Open output files for each phenotype
    std::vector<std::ofstream> out_files(n_pheno);
    for(int i = 0; i < n_pheno; ++i) {
        std::string fname = out_file + ".pheno" + std::to_string(i+1);
        out_files[i].open(fname, std::ios::out);
        // write header
        out_files[i] << "Chromosome\tPredictor\tBasepair\tA1\tA2\tBeta\tSE\tChisq\tP\tN\tMAF\tMiss\n";
    }

    for(int b = 0; b < nr_blocks; ++b){

        Rcout << "Performing linear regression for block " << b << "\n";

        // read in the genotype - need to optimize the matrix later
        int block_start = b * block_size;
        int block_end = std::min(n_snps - 1, (b + 1) * block_size - 1); // make sure you don't exceed n_snps
        int n_snps_block = block_end - block_start + 1;

        Rcpp::IntegerMatrix geno_block = readBedBlock(
            bed_file, 
            n_inds, n_snps, 
            0, n_inds - 1,              // full range of individuals
            block_start, block_end      // full range of SNPs
        );
        Eigen::MatrixXd genotypes = Rcpp::as<Eigen::MatrixXd>(geno_block);

        // Compute MAF and missingness
        Rcpp::NumericVector maf = computeMAF(geno_block);
        Rcpp::NumericVector miss = computeMissingness(geno_block);

        // Standardize columns - for one parameter alpha
        for (int col = 0; col < genotypes.cols(); ++col) {
            
            Eigen::VectorXd v = genotypes.col(col);
           
            double mean = v.mean();
            double stddev = std::sqrt((v.array() - mean).square().sum() / (v.size() - 1));
            
            if (stddev > 0) {
                genotypes.col(col) = (v.array() - mean) / stddev;
            } else {
                genotypes.col(col).setZero(); // trivial column (all same value), set to zero
            }
        }

        for(int i = 0; i < n_pheno; ++i){
                
            Eigen::VectorXd y = Eigen::VectorXd::Map(pheno.column(i).begin(), pheno.nrow());

            // Pre-allocate SNP outside
            std::vector<double> x_nonmiss(n_inds);
            std::vector<double> y_nonmiss(n_inds);

            // now loop over the SNPs
            for (int s = 0; s < n_snps_block; ++s) {
             
                int n_used = 0;

                // Fill pre-allocated vectors
                for (int ind = 0; ind < n_inds; ++ind) {
                    int g = geno_block(ind, s);
                    if (g != -1 && !NumericVector::is_na(y(ind))) {
                        x_nonmiss[n_used] = genotypes(ind, s);
                        y_nonmiss[n_used] = y(ind);
                        n_used++;
                    }
                }

                // Use Eigen::Map on pre-allocated vectors
                Eigen::Map<Eigen::VectorXd> xs2(x_nonmiss.data(), n_used);
                Eigen::Map<Eigen::VectorXd> ys2(y_nonmiss.data(), n_used);

                double xTx = xs2.dot(xs2);
                double beta = 0.0, se = NA_REAL, chisq = NA_REAL, pval = NA_REAL;

                if (xTx > 0) {
                    double xTy = xs2.dot(ys2);
                    beta = xTy / xTx;

                    double rss = (ys2 - xs2 * beta).squaredNorm();
                    double residual_var = rss / (n_used - 2);

                    se    = std::sqrt(residual_var / xTx);
                    chisq = (beta / se) * (beta / se);
                    pval  = 1.0 - R::pchisq(chisq, 1, 1, 0);   // upper tail
                }
                
                // Write output line
                out_files[i]
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