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
  K <- strata$K
  ids <- strata$y[, 1]

  # Compute eigenvector transformation
  trans <- eigen(gencov[1:K, 1:K])$vectors[, 1]

  # Compute medians for smoothing stratification variables
  group_table <- table(strata$info$groups)
  cum <- cumsum(group_table)
  medians <- (cum - 0.5 * (cum - c(0, cum[-length(cum)]))) / max(cum)

  # Compute medians for scaled stratification variable
  strat_scale <- as.numeric(scale(strata$info[, 3]))
  medians_obs <- unlist(lapply(1:K, function(x) mean(strat_scale[which(strata$info$groups == x)], na.rm = TRUE)))

  # Initialize transformed phenotype with controls
  trans_pheno <- cbind(ids, ids, 0)
  control_ids <- strata$y[strata$y[, 3] == 0, 1]

  if(strata$sparse){
    # For sparse case: use transformation weights on residualized case phenotypes
    for (k in 1:K) {
      group_ids <- strata$info[strata$info$groups == k, 1]
      trans_pheno[ids %in% group_ids, 3] <- trans[k]
    }

  } else {
    # Smooth medians through transformation
    #fit <- smooth.spline(medians, as.numeric(trans), spar = 0.2)
    #trans_pred <- predict(fit, strata$info$order)$y
    #names(trans_pred) <- strata$info[, 1]
    
    # Alternative, smooth directly through values
    fit <- smooth.spline(medians_obs, as.numeric(trans), spar = 0.2)
    trans_pred <- predict(fit, strat_scale)$y

    # Make the transformed phenotype (27-01-26)
    names(trans_pred) <- strata$info[, 1]

    # Scale phenotype to 0 / 1, then multiply with weight
    for (k in 1:K) {
      stratum <- strata[[k]]
      stratum[, 3] <- stratum[, 3] - mean(stratum[stratum[, 1] %in% control_ids, 3], na.rm = T)
      stratum[, 3] <- stratum[, 3] / mean(stratum[!stratum[, 1] %in% control_ids, 3], na.rm = T)

      weights <- trans_pred[match(stratum[, 1], names(trans_pred))]
      weights[is.na(weights)] <- mean(weights, na.rm = T)

      trans_pheno[match(stratum[, 1], trans_pheno[, 1]), 3] <- stratum[, 3] * weights
    }

    trans_pheno[trans_pheno[, 1] %in% control_ids, 3] <- trans_pheno[trans_pheno[, 1] %in% control_ids, 3] / K
    
    # Make the transformed phenotype (26-01-26)
    #names(trans_pred) <- strata$info[, 1]
    #idx <- match(names(trans_pred), trans_pheno[, 1])
    #valid <- !is.na(idx)
    #trans_pheno[idx[valid], 3] <- trans_pred[valid]
  } 

  # Add back cases with missing stratification variable
  if (length(strata$strat_miss) > 0) {
    # Compute mean phenotype value among cases (y == 1)
    case_ids <- strata$y[strata$y[, 3] == 1, 1]
    case_indices <- trans_pheno[, 1] %in% case_ids
    mean_case_pheno <- mean(as.numeric(trans_pheno[case_indices, 3]), na.rm = TRUE)

    trans_pheno_add <- cbind(strata$strat_miss, strata$strat_miss, mean_case_pheno)
    trans_pheno <- rbind(trans_pheno, trans_pheno_add)

    # Sort back to original order
    trans_pheno <- trans_pheno[match(strata$ids, trans_pheno[, 1]), ]
  }

  # Write to phenotype
  colnames(trans_pheno) = c("FID", "IID", "Pheno")
  write.table(trans_pheno, paste0(outfile, ".pheno"), quote = F, row = F)

  # Compute the inflation factor - need a2, gencor, and h2_Z
  h2_y <- gencov[K + 1, K + 1]
  h2_Z <- gencov[K + 2, K + 2]

  if(h2_y * h2_Z < 0){
    message(paste0("Can not compute expected inflation criterion due to negative h2 estimate of binary trait or stratification variable"))
  } else {
    # Prepare variables
    Z <- strata$Z[, 3]
    Z[is.na(Z)] <- mean(Z, na.rm = TRUE)
    
    y <- strata$y[, 3]
    y[is.na(y)] <- mean(y, na.rm = TRUE)
    
    # Scale variables
    Z <- as.numeric(scale(Z))
    y <- as.numeric(scale(y))
    Y_trans <- as.numeric(scale(trans_pheno[, 3]))
    
    # Fit regression to compute a2
    fit <- lm(Y_trans ~ y + Z - 1)
    coefs <- fit$coefficients
    a2 <- coefs["Z"]
    
    # Compute inflation criterion using genetic correlation
    rg <- gencov[K + 1, K + 2] / sqrt(h2_y * h2_Z)
    exp_inflation <- a2^2 * (1 - rg^2) * h2_Z
    
    message(sprintf("Expected inflation criterion is %.4f", exp_inflation))

    if (exp_inflation > 0.01) {
      warning("Inflation criterion is greater than 0.01.")
    }
  }

  invisible(trans_pheno)
}