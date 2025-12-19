// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <Rcpp.h>
#include "readBedBlock.h"
#include "geno_utils.h"
#include <vector>
#include <random>
#include <algorithm>

using namespace Rcpp;

using Eigen::Matrix2d;
using Eigen::Vector2d;

// [[Rcpp::export]]
Rcpp::NumericMatrix he_multi_part(const std::string& filename, const SEXP pheno_mat, int block_size) {
    
    // Get dimensions of .bim and .fam file, and read genotype ids
    int n_snps = count_lines(filename + ".bim");
    List fam = read_fam_file(filename);
    CharacterVector geno_iid = fam["iid"];
    int n_total_inds = geno_iid.size();

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

    // Output matrix
    Rcpp::NumericMatrix gencors(n_pheno, n_pheno);

    // Block size
    int nr_blocks = n_snps / block_size;

    // Set up random vectors
    std::mt19937 rng(12345);
    std::normal_distribution<double> norm(0,1);
    int nmcmc = 20;

    for(int b = 0; b < nr_blocks; ++b){

        if (b % 100 == 1) {
            Rcpp::Rcout << "Processing genotype block " << b << "/" << nr_blocks << "\n";
        }
    
        int block_start = b * block_size;
        int block_end = std::min(n_snps - 1, (b + 1) * block_size - 1);
        int n_snps_block = block_end - block_start + 1;

        // Read full genotype block
        IntegerMatrix geno_block_full = readBedBlock(
            filename + ".bed", 
            n_total_inds, n_snps, 
            0, n_total_inds - 1,              // full range of individuals
            block_start, block_end      // full range of SNPs
        );

        // Subset to matched individuals
        IntegerMatrix geno_block(n_inds, n_snps_block);
        for (int i = 0; i < n_inds; ++i) {
            for (int j = 0; j < n_snps_block; ++j) {
                geno_block(i, j) = geno_block_full(geno_keep[i], j);
            }
        }

        Eigen::MatrixXd X = Rcpp::as<Eigen::MatrixXd>(geno_block);

        // Compute MAF and weights - ignore for now
        // Rcpp::NumericVector maf = computeMAF(geno_block);
        // Rcpp::NumericVector weights(maf.size());

        // for (int j = 0; j < weights.size(); j++) {
        //    weights[j] = std::pow(maf[j] * (1 - maf[j]), 1 - 0.5);
        // }

        // Standardize columns - for one parameter alpha
        for (int col = 0; col < X.cols(); ++col) {
            
            Eigen::VectorXd v = X.col(col);
           
            // double w = weights[col]; 
            // v*=w;

            double mean = v.mean();
            double stddev = std::sqrt((v.array() - mean).square().sum() / (v.size() - 1));
            
            if (stddev > 0) {
                X.col(col) = (v.array() - mean) / stddev;
            } else {
                // trivial column (all same value), set to zero
                X.col(col).setZero();
            }
        }

        // Set traces for estimation
        double tr_K = n_inds;
        double tr_KK = 0.0;

        for (int r = 0; r < nmcmc; ++r) {
            Eigen::VectorXd z(n_inds);
            for (int i = 0; i < n_inds; ++i) z[i] = norm(rng);

            Eigen::VectorXd Kz = X * (X.transpose() * z);
            tr_KK += Kz.squaredNorm() / (nmcmc * n_snps_block * n_snps_block);
        }

        // Haseman-Elston equations
        for(int i = 0; i < n_pheno; ++i){
            for(int j = i; j < n_pheno; ++j){
                
                Eigen::VectorXd yi = Eigen::Map<Eigen::VectorXd>(pheno.column(i).begin(), n_inds);
                Eigen::VectorXd yj = Eigen::Map<Eigen::VectorXd>(pheno.column(j).begin(), n_inds);

                double yKyt = (X.transpose() * yi).dot(X.transpose() * yj) / n_snps_block;
                double yyt  = yi.dot(yj);

                Eigen::Matrix2d A;
                A << tr_KK, tr_K,
                     tr_K,  n_inds;

                Eigen::Vector2d B(yKyt, yyt);
                Eigen::Vector2d sol = A.colPivHouseholderQr().solve(B);

                gencors(i, j) += sol[0];
                if (i != j) gencors(j, i) += sol[0];
            }
        }
    }

    return gencors;
}