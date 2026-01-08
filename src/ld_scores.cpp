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
Eigen::MatrixXd convert_and_impute(const IntegerMatrix &geno, const std::vector<int> &inds_subset) {
  int N = inds_subset.size();
  int M = geno.ncol();
  Eigen::MatrixXd G(N, M);

  for(int c = 0; c < M; c++) {
      Eigen::ArrayXd col(N);
      double sum = 0.0;
      int cnt = 0;

      // One-pass imputation
      for(int r = 0; r < N; r++) {
        int g = geno(inds_subset[r], c);
          if(g >= 0) {
              col(r) = g;
              sum += g;
              cnt++;
          } else {
              col(r) = NAN;  // temporarily mark missing
          }
      }

      double mean = (cnt > 0) ? sum / cnt : 0.0;
      for(int r = 0; r < N; r++) {
          if(std::isnan(col(r))) col(r) = mean;
      }

      // Standardize
      double col_mean = col.mean();
      double col_sd = std::sqrt((col - col_mean).square().sum() / (N - 1));
      if(col_sd > 0) col = (col - col_mean) / col_sd;
      else col.setZero();

      G.col(c) = col.matrix();
  }

  return G;
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
DataFrame computeLDscoresFromBED(std::string file_prefix, IntegerVector geno_set) {

  // Read BIM and FAM files
  List bim = read_bim_file(file_prefix);
  List fam = read_fam_file(file_prefix);

  IntegerVector chr = bim["chr"];
  IntegerVector pos = bim["pos"];
  CharacterVector snp = bim["snp"];

  int n_snp = chr.size();
  int n_ind_total = as<CharacterVector>(fam["iid"]).size();

  // convert R indices to 0-based C++ indices
  std::vector<int> inds_subset(geno_set.size());
  for(int i = 0; i < geno_set.size(); ++i) inds_subset[i] = geno_set[i] - 1;

  int n_ind = inds_subset.size();
  Rcout << "Number of individuals = " << n_ind << ", number of SNPs = " << n_snp << "\n";

  // Determine block size
  int max_block = max_snps_in_1mb_block(chr, pos, 1e6);
  Rcout << "Largest 1Mb SNP block = " << max_block << "\n";
  int read_block_size = 2 * max_block;

  NumericVector ldscore(n_snp, 0.0);
  int start_snp = 0;
  static int last_printed = -1000; // Track last printed value

  // Strategy: read in blocks at a time to reduce memory. First find the maximum index distance that covers 1Mb? 
  while(start_snp < n_snp) {

    if (start_snp / 1000 > last_printed / 1000) {
      Rcout << "Reading SNP " << (start_snp / 1000) * 1000 << " / " << n_snp << "\n";
      last_printed = (start_snp / 1000) * 1000;
    }

    // Rcout << "Reading SNP " << start_snp << "/" << n_snp << "\n";

    // Determine the indices to compute LD for in this block
    int end_snp = std::min(start_snp + read_block_size - 1, n_snp - 1);
    int n_block = end_snp - start_snp + 1;

    int compute_start, compute_end;

    if(start_snp == 0 && end_snp == n_snp - 1) {
        // Single block covering all SNPs
        compute_start = 0;
        compute_end = end_snp;
    } else if(start_snp == 0) {
        // First block
        compute_start = 0;
        compute_end = std::min(start_snp + 3 * max_block / 2 - 1, n_snp - 1) - start_snp;
    } else if(end_snp == n_snp - 1) {
        // Last block
        compute_start = max_block / 2;
        compute_end = end_snp - start_snp;
    } else {
        // Middle block
        compute_start = max_block / 2;
        compute_end = std::min(3 * max_block / 2 - 1, end_snp - start_snp);
    }

    // Read BED block
    IntegerMatrix geno = readBedBlock(file_prefix + ".bed",
                                      n_ind_total, n_snp,
                                      0, n_ind_total - 1,
                                      start_snp, end_snp);

    Eigen::MatrixXd G = convert_and_impute(geno, inds_subset);

    // Compute correlation matrix for the block
    Eigen::MatrixXd cor = (G.transpose() * G) / (n_ind - 1);

    // Compute LD scores
    for(int j = compute_start; j <= compute_end; ++j) {
      int global_j = start_snp + j;
      double sum_r2 = 0.0;

      // Only consider SNPs on same chromosome within 1Mb
      for(int k = 0; k < n_block; ++k) {
          int global_k = start_snp + k;
          if(chr[global_j] != chr[global_k]) continue;
          if(std::abs(pos[global_j] - pos[global_k]) > 1e6) continue;

          double r = cor(j, k);
          sum_r2 += r * r;
      }

      ldscore[global_j] = sum_r2;
    }

    // Slide forward by max_block
    start_snp += max_block;
  }

  return DataFrame::create(
    Named("Predictor") = snp,
    Named("Tagging") = ldscore,
    _["stringsAsFactors"] = false
  );
}