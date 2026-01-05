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
  trans = eigen(gencov[1:K,1:K])$vectors[,1]

  # Compute medians for smoothing stratification variables
  cum = cumsum(table(strata$info$groups))
  medians = (cum - 1/2 * (cum - c(0, cum[-length(cum)])))/max(cum)

  # Smooth medians through transformation
  fit <- smooth.spline(medians, as.numeric(trans), spar = 0.2)
  trans_pred <- predict(fit, strata$info$order)$y
  names(trans_pred) = strata$info[,1]

  # Make the transformed phenotype
  trans_pheno = cbind(ids, ids, 0)

  idx = match(names(trans_pred), trans_pheno[,1])
  valid <- !is.na(idx)

  trans_pheno[idx[valid],3] = trans_pred[valid]
  colnames(trans_pheno) = c("FID", "IID", "Pheno")

  # Write to phenotype
  write.table(trans_pheno, paste0(outfile, ".pheno"), quote = F, row = F)

  # Also compute the inflation factor
  a2 = cor(strata$y[,3], trans_pheno[,3])
  exp_inflation = a2^2 * (1 - gencov[K+1, K+2]^2) * gencov[K+2, K+2]
  message(paste0("Expected inflation criterion is ",round(exp_inflation, 4)))
  if(exp_inflation > 0.01) warning("Inflation criterion is greater than 0.01.")
}