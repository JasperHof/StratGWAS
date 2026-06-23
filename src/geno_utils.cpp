#include <Rcpp.h>
#include "geno_utils.h"
#include "readBedBlock.h"
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
using namespace Rcpp;

// Function for computing MAF
NumericVector computeMAF(const IntegerMatrix &geno) {
    int n_snp = geno.ncol();
    int n_ind = geno.nrow();

    NumericVector maf(n_snp);

    for(int j = 0; j < n_snp; j++) {
        int sum = 0, count = 0;
        for(int i = 0; i < n_ind; i++) {
            int g = geno(i,j);
            if(g >= 0) {       // skip missing (-1)
                sum += g;
                count++;
            }
        }
        double p = sum / (2.0 * count);
        maf[j] = std::min(p, 1-p);
    }
    
    return maf;
}

// Function for computing missingness (coded as -1)
NumericVector computeMissingness(const IntegerMatrix &geno) {
    int n_snp = geno.ncol();
    int n_ind = geno.nrow();

    NumericVector miss(n_snp);

    for(int j = 0; j < n_snp; j++) {
        int missing = 0;
        for(int i = 0; i < n_ind; i++) {
            if(geno(i,j) == -1)
                missing++;
        }
        miss[j] = missing / double(n_ind);
    }
    return miss;
}

// [[Rcpp::export]]
List read_bim_file(const std::string &filename) {
    // Construct the BIM filename
    std::string bim_file = filename + ".bim";

    // Call R's read.table() function
    Function read_table("read.table");
    List bim_data = read_table(
        Named("file") = bim_file,
        Named("header") = false,
        Named("stringsAsFactors") = false,
        Named("colClasses") = CharacterVector::create("integer", "character", "numeric", "integer", "character", "character")
    );

    // Extract columns
    CharacterVector snp = bim_data[1];
    IntegerVector chr = bim_data[0];
    IntegerVector pos = bim_data[3];
    CharacterVector a1 = bim_data[4];
    CharacterVector a2 = bim_data[5];

    // Return as a list (or process further in C++)
    return List::create(
        Named("snp") = snp,
        Named("chr") = chr,
        Named("pos") = pos,
        Named("a1") = a1,
        Named("a2") = a2
    );
}

// [[Rcpp::export]]
Rcpp::List read_fam_file(const std::string &filename) {

    std::string fam_file = filename + ".fam";

    Function read_table("read.table");
    List fam_data = read_table(
        Named("file") = fam_file,
        Named("header") = false,
        Named("stringsAsFactors") = false,
        Named("colClasses") = CharacterVector::create(
            "character",  // FID
            "character",  // IID
            "character",  // father
            "character",  // mother
            "integer",    // sex
            "numeric"     // phenotype
        )
    );

    return List::create(
        Named("fid")   = fam_data[0],
        Named("iid")   = fam_data[1],
        Named("father")= fam_data[2],
        Named("mother")= fam_data[3],
        Named("sex")   = fam_data[4],
        Named("pheno") = fam_data[5]
    );
}

// Function just to read dimension of .bim and .fam files
int count_lines(const std::string& filename) {
    std::ifstream file(filename.c_str());
    if (!file.is_open()) {
        Rcpp::stop("Cannot open file: " + filename);
    }

    int count = 0;
    std::string line;
    while (std::getline(file, line)) {
        ++count;
    }
    return count;
}

// [[Rcpp::export]]
Rcpp::NumericVector compute_maf_all(const std::string& filename) {
 
    List bim = read_bim_file(filename);
    List fam = read_fam_file(filename);
 
    int n_snps       = Rcpp::as<Rcpp::CharacterVector>(bim["snp"]).size();
    int n_total_inds = Rcpp::as<Rcpp::CharacterVector>(fam["iid"]).size();
 
    NumericVector maf_all(n_snps);
 
    int block_size = 256;
    int n_blocks   = (n_snps + block_size - 1) / block_size;
 
    for (int b = 0; b < n_blocks; ++b) {
 
        int snp_start  = b * block_size;
        int snp_end    = std::min(n_snps - 1, snp_start + block_size - 1);
        int n_snps_blk = snp_end - snp_start + 1;
 
        IntegerMatrix geno_block = readBedBlock(
            filename + ".bed",
            n_total_inds, n_snps,
            0, n_total_inds - 1,
            snp_start, snp_end
        );
 
        NumericVector maf_blk = computeMAF(geno_block);
 
        for (int s = 0; s < n_snps_blk; ++s)
            maf_all[snp_start + s] = maf_blk[s];
    }
 
    return maf_all;
}