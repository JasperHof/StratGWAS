#' Update eigendecomposition
#'
#' Update eigendecomposition to meet the inflation criterion
#'
#' @useDynLib StratGWAS, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' 
#' @param cor_g Matrix containing genetic correlations
#' @param cor_e Matrix containing environmental correlations
#' @param multi Input design matrix
#' @return Returns covariance matrix of the strata
#' @export
find_optimal_scaling <- function(cor_g, cor_e, multi, max_criterion = 0.01, 
                                 tolerance = 0.001, max_iter = 20, verbose = TRUE) {
  
  # Function for multiplying auxiliary variable heritabilities by (0 < factor <= 1)
  downscale_auxiliary_heritability <- function(cor_g, scaling_factor = 0.5) {
    cor_g_scaled <- cor_g
    n <- nrow(cor_g)
    
    if(n < 2) return(cor_g)
    
    # Keep first variable (binary trait) unchanged
    for(i in 2:n) {
      original_h2 <- cor_g[i, i]
      scaled_h2 <- original_h2 * scaling_factor
      
      # Scale the correlations proportionally
      # If h2_new = alpha * h2_old, then Cov_new = sqrt(alpha) * Cov_old
      scale_sqrt <- sqrt(scaling_factor)
      cor_g_scaled[i, ] <- cor_g_scaled[i, ] * scale_sqrt
      cor_g_scaled[, i] <- cor_g_scaled[, i] * scale_sqrt
      cor_g_scaled[i, i] <- scaled_h2  # Restore diagonal
    }
    
    return(cor_g_scaled)
  }
  
  # Function to compute criterion for a given scaling factor
  compute_criterion <- function(scaling) {
    cor_g_scaled <- downscale_auxiliary_heritability(cor_g, scaling)
    
    # Perform eigendecomposition
    P <- eigen(cor_e)
    S_inv <- P$vectors %*% diag(1/sqrt(P$values)) %*% t(P$vectors)
    M <- S_inv %*% as.matrix(cor_g_scaled) %*% S_inv
    V <- eigen(M)$vectors
    U1 <- t(V) %*% S_inv
    weights <- U1[1,]
    
    # Compute transformed phenotype
    y_trans <- as.numeric(scale(multi %*% weights))
    y <- as.numeric(scale(multi[, 1]))
    z <- as.numeric(scale(multi[, 2]))
    
    # Fit regression
    fit <- lm(y_trans ~ y + z - 1)
    coefs <- fit$coefficients
    a2 <- coefs["z"]
    
    # Compute genetic correlation
    rg <- cor_g_scaled
    hers <- diag(cor_g_scaled)
    
    for(k in 1:nrow(rg)) {
      if(hers[k] > 0) {
        rg[k, ] <- rg[k, ] / sqrt(hers[k])
        rg[, k] <- rg[, k] / sqrt(hers[k])
      }
    }
    
    gamma <- rg[1, 2]
    h2_Z <- cor_g_scaled[2, 2]
    
    if(is.na(gamma) || is.na(a2)) {
      criterion <- NA
    } else {
      criterion <- a2^2 * (1 - gamma^2) * h2_Z
    }
    
    return(list(criterion = criterion, scaling = scaling, weights = weights))
  }
  
  # NEW: First check if scaling is needed (with scaling_factor = 1.0)
  if(verbose) message("Checking if inflation constraint is needed...")
  
  initial_result <- compute_criterion(scaling = 1.0)
  
  if(!is.na(initial_result$criterion)) {
    if(verbose) {
      message(sprintf("Current criterion (unscaled): %.6f", initial_result$criterion))
    }
    
    if(initial_result$criterion <= max_criterion) {
      if(verbose) {
        message("Criterion already satisfies constraint. No scaling needed.")
      }
      # Return unscaled result
      return(list(
        optimal_scaling = 1.0,
        criterion = initial_result$criterion,
        weights = initial_result$weights,
        cor_g_scaled = cor_g,
        scaling_applied = FALSE
      ))
    } else {
      if(verbose) {
        message("Criterion exceeds threshold. Applying optimal scaling...")
      }
    }
  } else {
    if(verbose) {
      message("Cannot compute initial criterion (NA values). Proceeding with scaling...")
    }
  }
  
  # Binary search to find optimal scaling
  lower <- 0.001
  upper <- 1.0
  
  for(iter in 1:max_iter) {
    mid <- (lower + upper) / 2
    result <- compute_criterion(mid)
    
    if(is.na(result$criterion)) {
      if(verbose) cat("Criterion is NA, increasing scaling\n")
      lower <- mid
      next
    }
    
    if(verbose) {
      cat(sprintf("Iter %d: scaling=%.4f, criterion=%.6f\n", 
                  iter, mid, result$criterion))
    }
    
    if(abs(result$criterion - max_criterion) < tolerance) {
      if(verbose) cat("Converged!\n")
      return(list(
        optimal_scaling = mid,
        criterion = result$criterion,
        weights = result$weights,
        cor_g_scaled = downscale_auxiliary_heritability(cor_g, mid),
        scaling_applied = TRUE
      ))
    }
    
    if(result$criterion > max_criterion) {
      # Need to decrease scaling (reduce auxiliary variance)
      upper <- mid
    } else {
      # Can increase scaling
      lower <- mid
    }
  }
  
  if(verbose) warning("Did not converge within max iterations")
  final_result <- compute_criterion(mid)
  return(list(
    optimal_scaling = mid,
    criterion = final_result$criterion,
    weights = final_result$weights,
    cor_g_scaled = downscale_auxiliary_heritability(cor_g, mid),
    scaling_applied = TRUE
  ))
}