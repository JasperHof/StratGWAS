#include <RcppEigen.h>
#include "geno_utils.h"
#include "readBedBlock.h"
#include <vector>
#include <fstream>
#include <cmath>
#include <string>
#include <bitset>

using namespace Rcpp;

// Helper function that determines how large the block size should be
int max_snps_in_1mb_block(IntegerVector chr, IntegerVector pos, int window_size = 1e6) {
  int n = chr.size();
  int min_diff = 0;
  int i = 1;

  while(i < n - 1 && min_diff < 1e6){
    IntegerVector diffs = pos[Range(1 + i, n - 1)] - pos[Range(0, n - 2 - i)];
    LogicalVector same_chr = chr[Range(1 + i, n - 1)] == chr[Range(0, n - 2 - i)];

    IntegerVector diffs_same_chr = diffs[same_chr];

    if (diffs_same_chr.size() == 0)
        break;

    min_diff = min(diffs_same_chr);
    ++i;
  }

  return i;
}

// Helper function to impute missing values and standardise
Eigen::MatrixXd convert_and_impute(const IntegerMatrix &geno) {
  int N = geno.nrow();
  int M = geno.ncol();
  Eigen::MatrixXd G(N, M);

  for (int c = 0; c < M; c++) {
      // impute mean genotype
      double sum = 0;
      int cnt = 0;
      for (int r = 0; r < N; r++) {  
          int g = geno(r, c);
          if (g >= 0) { sum += g; cnt++; }
      }
      double mean = (cnt > 0 ? sum / cnt : 0.0);
      for (int r = 0; r < N; r++) {
          int g = geno(r, c);
          G(r, c) = (g >= 0 ? g : mean);
      }

      // standardize
      Eigen::VectorXd v = G.col(c);
      double sd = std::sqrt((v.array() - v.mean()).square().sum() / (N - 1));
      if (sd > 0) G.col(c) = (v.array() - v.mean()) / sd;
      else G.col(c).setZero();
  }

  return G;
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
NumericVector computeLDscoresFromBED(std::string file_prefix, int n_ind, int n_snp) {
  
  // First load BIM file
  Rcout << "Reading BIM file " << "\n";

  // Call the new function and extract vectors
  List bim = read_bim_file(file_prefix);
  CharacterVector snp = bim["snp"];
  IntegerVector chr = bim["chr"];
  IntegerVector pos = bim["pos"];
  CharacterVector a1 = bim["a1"];
  CharacterVector a2 = bim["a2"];

  Rcout << "Read BIM file " << "\n";

  // Compute maximum number of SNPs within 1Mb to determine sliding block size
  int max_block = max_snps_in_1mb_block(chr, pos, 1e6);
  Rcout << "Largest 1Mb SNP block = " << max_block << "\n";

  int read_block_size = 2 * max_block;
  NumericVector ldscore(n_snp);

  int start_snp = 0;

  // Strategy: read in blocks at a time to reduce memory. First find the maximum index distance that covers 1Mb? 
  while(start_snp < n_snp) {
    int end_snp = std::min(start_snp + read_block_size - 1, n_snp - 1);

    // Determine the indices to compute LD for in this block
    int compute_start, compute_end;

    if(start_snp == 0 && end_snp == n_snp - 1) {
        // Single block covering all SNPs
        compute_start = 0;
        compute_end = end_snp;
    } else if(start_snp == 0) {
        // first block
        compute_start = 0;
        compute_end = std::min(start_snp + 3 * max_block / 2 - 1, n_snp - 1) - start_snp;
    } else if(end_snp == n_snp - 1) {
        // last block
        compute_start = max_block / 2;
        compute_end = end_snp - start_snp;
    } else {
        // middle block
        compute_start = max_block / 2;
        compute_end = std::min(3 * max_block / 2 - 1, end_snp - start_snp);
    }

    // Read BED block
    IntegerMatrix geno = readBedBlock(file_prefix + ".bed",
                                      n_ind, n_snp,
                                      0, n_ind - 1,
                                      start_snp, end_snp);

    Eigen::MatrixXd G = convert_and_impute(geno);

    // Compute LD scores for this block
    for(int j = compute_start; j <= compute_end; ++j){
        int global_j = start_snp + j;
        double sum_r2 = 0.0;

        for(int k = 0; k < G.cols(); ++k){
            int global_k = start_snp + k;
            if(chr[global_j] != chr[global_k]) continue;
            if(std::abs(pos[global_j] - pos[global_k]) > 1e6) continue;

            double r = G.col(j).dot(G.col(k)) / (n_ind - 1);
            sum_r2 += r * r;
        }
        ldscore[global_j] = sum_r2;
    }

    // Slide forward by max_block
    start_snp += max_block;
  }

  return ldscore;
}