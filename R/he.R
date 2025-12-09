#' Haseman-Elston regression
#'
#' A short description of what the function does.
#'
#' @param genotypes Numeric matrix of genotypes (individuals x SNPs)
#' @param phenotype Numeric vector of phenotypes
#' @return Numeric vector of length 2: h2 and e2
#' @export
he <- function(genotypes, phenotype) {
  # input checks
  check_genotype(genotypes)
  check_pheno(phenotype, nrow(genotypes))
  
  # call C++ backend
  he(genotypes, phenotype)
}