#' LD score regression
#'
#' Compute genetic covariance of subgroups using LD score regression
#'
#' @useDynLib StratGWAS, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' 
#' @param strata An object returned from stratify()
#' @param gencov Genetic covariance matrix
#' @param outfile Name of output file for phenotype
#' @return Returns covariance matrix of the strata
#' @export
transform <- function(strata, gencov, outfile) {
  
  # <performs some checks here> #

  # Read in data as multivariate phenotype
  K = strata$K
  ids = strata$y[,1]
  trans = eigen(gencov)$vectors[,1]

  # Compute medians for smoothing stratification variables
  cum = cumsum(table(strata$info$groups))
  medians = (cum - 1/2 * (cum - c(0, cum[-length(cum)])))/max(cum)

  # Smooth medians through transformation
  fit <- smooth.spline(medians, as.numeric(trans), spar = 0.2)
  trans_pred <- predict(fit, strata$info$order)$y
  names(trans_pred) = strata$info[,1]

  trans_pheno = cbind(ids, ids, 0)
  trans_pheno[match(names(trans_pred), trans_pheno[,1]),3] = trans_pred
  colnames(trans_pheno) = c("FID", "IID", "Pheno")

  # Write to phenotype
  write.table(trans_pheno, paste0(outfile, ".pheno"), quote = F, row = F)
}