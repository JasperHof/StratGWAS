#ifndef READBEDBLOCK_H
#define READBEDBLOCK_H

#include <Rcpp.h>
#include <string>

// [[Rcpp::depends(Rcpp)]]

// Declaration of the function
Rcpp::IntegerMatrix readBedBlock(std::string filename, int n_ind, int n_snp, 
                                 int start_ind, int end_ind, int start_snp, int end_snp);

#endif // READBEDBLOCK_H
