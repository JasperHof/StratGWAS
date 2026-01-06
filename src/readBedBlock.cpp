#include <Rcpp.h>
#include <vector>
#include <fstream>
#include <cmath>
#include <string>
#include <bitset>

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::IntegerMatrix readBedBlock(std::string filename, int n_ind, int n_snp, 
                               int start_ind, int end_ind, int start_snp, int end_snp) {
  
  std::ifstream bedFile(filename, std::ios::binary);
  if(!bedFile.is_open()) {
    stop("Error opening file " + filename);
  }
  
  // Read magic numbers
  unsigned char byte1, byte2, byte3;
  bedFile.read(reinterpret_cast<char*>(&byte1), 1);
  bedFile.read(reinterpret_cast<char*>(&byte2), 1);
  bedFile.read(reinterpret_cast<char*>(&byte3), 1);
  
  if(byte1 != 108 || byte2 != 27 || byte3 != 1) {
    stop("BED file header invalid or unsupported format (only SNP-major).");
  }
  
  int bytes_per_snp = static_cast<int>(std::ceil(n_ind / 4.0));
  std::vector<unsigned char> buffer(bytes_per_snp);
  
  int n_snps_block = end_snp - start_snp + 1;
  int n_inds_block = end_ind - start_ind + 1;
  
  IntegerMatrix genoMatrix(n_inds_block, n_snps_block);
  
  // Seek once to the start of the SNP block
  std::streampos start_pos = 3 + static_cast<std::streampos>(start_snp) * bytes_per_snp;
  bedFile.seekg(start_pos, std::ios::beg);
  
  if (!bedFile) {
      stop("Error seeking to SNP position " + std::to_string(start_snp));
  }

  for(int snp = start_snp; snp <= end_snp; ++snp) {

    // Read one SNP worth of data
    bedFile.read(reinterpret_cast<char*>(buffer.data()), bytes_per_snp);
        
    if (!bedFile) {
        stop("Error reading SNP " + std::to_string(snp));
    }
    
    // Decode genotypes for the requested individuals
    for (int ind = start_ind; ind <= end_ind; ++ind) {
        int byte_index = ind / 4;
        int shift = (ind % 4) * 2;
        
        unsigned char geno = (buffer[byte_index] >> shift) & 0b11;
        
        int g;
        switch (geno) {
            case 0b00: g = 2; break;  // Homozygous alternative
            case 0b10: g = 1; break;  // Heterozygous
            case 0b11: g = 0; break;  // Homozygous reference
            case 0b01: g = -1; break; // Missing
            default: g = -2; break;   // Should never happen
        }
        
        genoMatrix(ind - start_ind, snp - start_snp) = g;
    }

    /*
    bedFile.seekg(3 + snp * bytes_per_snp, std::ios::beg);
    bedFile.read(reinterpret_cast<char*>(buffer.data()), bytes_per_snp);
    
    for(int ind = start_ind; ind <= end_ind; ++ind) {
      int byte_index = ind / 4;
      int shift = (ind % 4) * 2;
      unsigned char geno = (buffer[byte_index] >> shift) & 0b11;
      
      int g;
      switch(geno) {
        case 0b00: g = 2; break;
        case 0b10: g = 1; break;  
        case 0b11: g = 0; break;
        case 0b01: g = -1; break;  // missing
        default: g = -2;
      }
      genoMatrix(ind - start_ind, snp - start_snp) = g;
    }
    */
  }

  bedFile.close();
  return genoMatrix;
}
