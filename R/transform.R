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

  if(strata$sparse){    
    # Make the transformed phenotype - directly from eigendecomposition
    trans_pheno = cbind(ids, ids, 0)
    for(k in 1:K) trans_pheno[ids %in% strata$info[strata$info$groups == k,1],3] <- trans[k]
  } else {
    # Smooth medians through transformation
    fit <- smooth.spline(medians, as.numeric(trans), spar = 0.2)
    trans_pred <- predict(fit, strata$info$order)$y
    names(trans_pred) = strata$info[,1]

    # Make the transformed phenotype
    trans_pheno = cbind(ids, ids, 0)
    idx = match(names(trans_pred), trans_pheno[,1])
    valid <- !is.na(idx)

    trans_pheno[idx[valid],3] = trans_pred[valid]
  } 

  # Add back cases with missing stratification variable
  if(length(strata$strat_miss) > 0){
    trans_pheno_add <- cbind(strata$strat_miss, strata$strat_miss,
                            mean(trans_pred[valid]))
    trans_pheno <- rbind(trans_pheno, trans_pheno_add)

    # Sort back to original order
    trans_pheno <- trans_pheno[match(strata$ids, trans_pheno[, 1]), ]
  }

  # Write to phenotype
  colnames(trans_pheno) = c("FID", "IID", "Pheno")
  write.table(trans_pheno, paste0(outfile, ".pheno"), quote = F, row = F)

  # Also compute the inflation factor - need a2, gencor, and h2_Z
  if(gencov[K+1, K+1] * gencov[K+2, K+2] < 0){
    message(paste0("Can not compute expected inflation criterion due to negative h2 estimate of binary trait or stratification variable"))
  } else {

    # Prepare variables
    Z = strata$Z[,3]
    Z[is.na(Z)] = mean(Z, na.rm=T)
    y = strata$y[,3]
    y[is.na(y)] = mean(y, na.rm=T)

    # Scale variables
    Z = as.numeric(scale(Z))
    y = as.numeric(scale(y))
    Y_trans = as.numeric(scale(trans_pheno[,3]))

    # Fit regression to compute a2
    fit = lm(Y_trans ~ y + Z - 1)
    coefs = fit$coefficients
    a2 = coefs["Z"]

    # Compute inflation criterion using genetic correlation
    rg = gencov[K+1, K+2] / sqrt(gencov[K+1, K+1] * gencov[K+2, K+2])
    exp_inflation = a2^2 * (1 - rg^2) * gencov[K+2, K+2]
    message(paste0("Expected inflation criterion is ",round(exp_inflation, 4)))
    if(exp_inflation > 0.01) warning("Inflation criterion is greater than 0.01.")
  }
}