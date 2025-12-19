#ifndef GENO_UTILS_H
#define GENO_UTILS_H

#include <Rcpp.h>

Rcpp::NumericVector computeMAF(const Rcpp::IntegerMatrix &geno);
Rcpp::NumericVector computeMissingness(const Rcpp::IntegerMatrix &geno);
Rcpp::List read_bim_file(const std::string &filename);
Rcpp::List read_fam_file(const std::string &filename);
int count_lines(const std::string& filename);

#endif