#' Implementation of SumHer for estimating SNP Heritability
#'
#' @param ss Summary statistics data.frame with columns: SNP, N, Chisq
#' @param ldscores Numeric vector of LD scores
#' @param fit_intercept Logical, whether to fit intercept (default TRUE)
#' @param tol Convergence tolerance (default 1e-6)
#' @param max_iter Maximum iterations (default 100)
#' @param weights Optional weights (default: 1/ldscores)
#' @return List with heritability estimates and likelihood
#' @export
sumher <- function(ss, ldscores, 
                    fit_intercept = TRUE,
                    tol = 1e-6,
                    max_iter = 100,
                    alpha = -0.25) {

  M <- length(ldscores)
  N_mean <- mean(ss$N)
  N_scaled <- ss$N / N_mean
  chisq_obs <- ss$Chisq
  MAF <- ss$MAF
  q <- (MAF * (1 - MAF))^(1 + alpha)
  q <- q / sum(q) 

  n_params <- 1 + fit_intercept
  theta <- rep(0, n_params)
  
  # Initial weights: D_inv = ldscores (approximately sum of r²)
  D_inv_0 <- pmax(ldscores, 1)  # Keep this fixed for likelihood comparison
  
  # Current weights (will be updated)
  D_inv <- D_inv_0
  
  # Starting expectations
  mu <- rep(1, M)

  # Compute log-likelihood using FIXED initial weights
  compute_loglik <- function(mu, chisq_obs, D_inv_fixed) {
    D <- 1 / D_inv_fixed
    d <- sum(D)
    lambda <- sum(D * (chisq_obs - mu)^2)
    loglik <- -d/2 * (log(2 * pi * lambda / d) + 1)
    return(loglik)
  }
  
  # Compute expectations from parameters
  compute_expectations <- function(theta, N_scaled, ldscores, q, M, fit_intercept) {
    mu <- rep(1, length(ldscores)) + theta[1] * M * N_scaled * q * ldscores
    if (fit_intercept) {
      mu <- mu + theta[2] * N_scaled
    }
    mu <- pmax(mu, 1e-6)
    return(mu)
  }
  

  # Initial log-likelihood
  loglik_old <- compute_loglik(mu, chisq_obs, D_inv_0)
  
  converged <- FALSE
  iter <- 0
  
  while (iter < max_iter && !converged) {
    iter <- iter + 1
    
    # Build design matrix weighted by CURRENT iteration weights
    W <- 1 / D_inv  # Current weights
    W_sqrt <- sqrt(W)
    
    X <- matrix(0, nrow = M, ncol = n_params)
    X[, 1] <- W_sqrt * M * N_scaled * q * ldscores
    if (fit_intercept) {
      X[, 2] <- W_sqrt * N_scaled
    }
    
    # Weighted response: sqrt(W) * (S - 1)
    y <- W_sqrt * (chisq_obs - 1)
    
    # Weighted least squares: solve X'X theta = X'y
    XtX <- crossprod(X)
    Xty <- crossprod(X, y)
    
    theta_new <- tryCatch({
      solve(XtX, Xty)
    }, error = function(e) {
      warning("Singular matrix, using previous theta")
      theta
    })
    
    # Update expectations
    mu <- compute_expectations(theta_new, N_scaled, ldscores, q, M, fit_intercept)
    
    # Update weights for NEXT iteration (but keep D_inv_0 for likelihood)
    # D_inv = ldscores * mu (where mu = expected value)
    D_inv <- pmax(ldscores, 1) * mu
    
    # Compute likelihood using FIXED initial weights
    loglik_new <- compute_loglik(mu, chisq_obs, D_inv_0)
    
    diff <- loglik_new - loglik_old
    theta <- theta_new
    loglik_old <- loglik_new
    
    h2_snp <- theta[1]
    
    if (abs(diff) < tol) {
      converged <- TRUE
    }
  }

  # Final estimates
  h2_snp <- theta[1]
  intercept <- 1
  if (fit_intercept) {
    intercept <- 1 + theta[2]
  }
  
  # Standard errors via inverse of X'X (at final weights)
  W <- 1 / D_inv
  W_sqrt <- sqrt(W)
  X <- matrix(0, nrow = M, ncol = n_params)
  X[, 1] <- W_sqrt * M * N_scaled * q * ldscores
  if (fit_intercept) {
    X[, 2] <- W_sqrt * N_scaled
  }
  
  vcov <- tryCatch({
    solve(crossprod(X))
  }, error = function(e) {
    matrix(NA, nrow = n_params, ncol = n_params)
  })
  
  se_h2 <- sqrt(vcov[1, 1])
  
  return(list(
    h2_snp = h2_snp,
    se_h2 = se_h2,
    intercept = intercept,
    loglik = loglik_old,
    theta = theta,
    converged = converged,
    iterations = iter,
    mu = mu
  ))
}


#' Implementation of SumHer for estimating genetic covariance
#'
#' @param ss1 Summary statistics for trait 1: data.frame with SNP, N, Z
#' @param ss2 Summary statistics for trait 2: data.frame with SNP, N, Z
#' @param ldscores Numeric vector of LD scores
#' @param h2_1 Heritability estimate for trait 1 (optional, will estimate if NULL)
#' @param h2_2 Heritability estimate for trait 2 (optional, will estimate if NULL)
#' @param sample_overlap Proportion of sample overlap (default 0)
#' @param fit_intercept Logical, fit intercept for overlap (default TRUE)
#' @param tol Convergence tolerance (default 1e-6)
#' @param max_iter Maximum iterations (default 100)
#' @return List with genetic correlation and covariance estimates
#' @export
sumher_cov <- function(ss1, ss2, ldscores,
                        h2_1 = NULL, h2_2 = NULL,
                        sample_overlap = 0,
                        fit_intercept = TRUE,
                        tol = 1e-6,
                        max_iter = 100,
                        alpha = -0.25) {
  
  # Input validation
  if (nrow(ss1) != nrow(ss2) || nrow(ss1) != length(ldscores)) {
    stop("Summary statistics and LD scores must match in length")
  }
  
  M <- length(ldscores)
  MAF <- ss1$MAF
  q <- (MAF * (1 - MAF))^(1 + alpha)
  q <- q / sum(q) 

  # Scale sample sizes
  N_mean_1 <- mean(ss1$N)
  N_mean_2 <- mean(ss2$N)
  N_scaled_1 <- ss1$N / N_mean_1
  N_scaled_2 <- ss2$N / N_mean_2
  
  # Product of Z-statistics (signed test statistics)
  Z_product <- (ss1$Beta / ss1$SE) * (ss2$Beta / ss2$SE)
  
  # Cross-product sample size (geometric mean)
  N_cross <- sqrt(ss1$N) * sqrt(ss2$N)
  N_cross_mean <- sqrt(N_mean_1) * sqrt(N_mean_2)
  N_cross_scaled <- N_cross / N_cross_mean
  
  # Estimate h2 for each trait if not provided
  if (is.null(h2_1)) {
    # cat("Estimating h2 for Trait 1...\n")
    res1 <- sumher(
      data.frame(SNP = ss1$Predictor, N = ss1$N, Chisq = ss1$Chisq, MAF = ss$MAF),
      ldscores = ldscores,
      fit_intercept = TRUE,
      tol = tol,
      max_iter = max_iter,
      alpha = alpha
    )
    h2_1 <- res1$h2_snp
  }
  
  if (is.null(h2_2)) {
    #cat("Estimating h2 for Trait 2...\n")
    res2 <- sumher(
      data.frame(SNP = ss2$Predictor, N = ss2$N, Chisq = ss2$Chisq, MAF = ss$MAF),
      ldscores = ldscores,
      fit_intercept = TRUE,
      tol = tol,
      max_iter = max_iter,
      alpha = alpha
    )
    h2_2 <- res2$h2_snp
  }
  
  # Initialize parameters
  n_params <- 1 + fit_intercept
  theta <- rep(0, n_params)
  
  # Initial weights
  D_inv_0 <- pmax(ldscores, 1)
  D_inv <- D_inv_0
  
  # Starting expectations (null model: no genetic covariance)
  # E[Z_A * Z_B] = c_AB (phenotypic correlation from sample overlap)
  mu <- rep(sample_overlap, M)
  
  # Compute log-likelihood
  compute_loglik <- function(mu, Z_product, D_inv_fixed) {
    D <- 1 / D_inv_fixed
    d <- sum(D)
    lambda <- sum(D * (Z_product - mu)^2)
    loglik <- -d/2 * (log(2 * pi * lambda / d) + 1)
    return(loglik)
  }
  
  # Compute expectations from parameters
  compute_expectations <- function(theta, N_cross_scaled, ldscores, q, M,
                                  sample_overlap, fit_intercept) {
    # E[Z_A * Z_B] = c_AB + h2_AB * N_cross * ldscores
    mu <- rep(sample_overlap, length(ldscores)) + 
          theta[1] * M * N_cross_scaled * q * ldscores
    
    if (fit_intercept) {
      # Intercept captures additional sample overlap effects
      mu <- mu + theta[2] * N_cross_scaled
    }
    
    return(mu)
  }
  
  # Initial log-likelihood
  loglik_old <- compute_loglik(mu, Z_product, D_inv_0)
  
  converged <- FALSE
  iter <- 0
  
  while (iter < max_iter && !converged) {
    iter <- iter + 1
    
    # Build design matrix weighted by CURRENT iteration weights
    W <- 1 / D_inv
    W_sqrt <- sqrt(W)
    
    X <- matrix(0, nrow = M, ncol = n_params)
    X[, 1] <- W_sqrt * M * N_cross_scaled * q * ldscores
    if (fit_intercept) {
      X[, 2] <- W_sqrt * N_cross_scaled
    }
    
    # Weighted response: sqrt(W) * (Z_A * Z_B - c_AB)
    y <- W_sqrt * (Z_product - sample_overlap)
    
    # Weighted least squares
    XtX <- crossprod(X)
    Xty <- crossprod(X, y)
    
    theta_new <- tryCatch({
      solve(XtX, Xty)
    }, error = function(e) {
      warning("Singular matrix, using previous theta")
      theta
    })
    
    # Update expectations
    mu <- compute_expectations(theta_new, N_cross_scaled, ldscores, q, M,
                               sample_overlap, fit_intercept)
    
    # Update weights for NEXT iteration
    # For cross-trait, weights depend on sqrt of individual trait expectations
    # Approximate as: D_inv ≈ ldscores * sqrt(mu_1 * mu_2)
    # Simpler approximation: use geometric mean of ldscores and absolute mu
    D_inv <- pmax(ldscores, 1) * pmax(abs(mu), 0.1)
    
    # Compute likelihood using FIXED initial weights
    loglik_new <- compute_loglik(mu, Z_product, D_inv_0)
    
    diff <- loglik_new - loglik_old
    theta <- theta_new
    loglik_old <- loglik_new
    
    # Convert to genetic covariance and correlation
    h2_AB <- theta[1]
    
    if (abs(diff) < tol) {
      converged <- TRUE
    }
  }
  
  # Final estimates
  h2_AB <- theta[1]
  rg <- NA_real_
  if (!is.na(h2_1) && !is.na(h2_2) && h2_1 > 0 && h2_2 > 0) {
    rg <- h2_AB / sqrt(h2_1 * h2_2)
  }
  
  intercept <- sample_overlap
  if (fit_intercept) {
    intercept <- sample_overlap + theta[2]
  }
  
  # Standard errors via inverse of X'X (at final weights)
  W <- 1 / D_inv
  W_sqrt <- sqrt(W)
  X <- matrix(0, nrow = M, ncol = n_params)
  X[, 1] <- W_sqrt * N_cross_scaled * ldscores
  if (fit_intercept) {
    X[, 2] <- W_sqrt * N_cross_scaled
  }
  
  vcov <- tryCatch({
    solve(crossprod(X))
  }, error = function(e) {
    matrix(NA, nrow = n_params, ncol = n_params)
  })
  
  # SE for h2_AB
  se_h2_AB <- sqrt(vcov[1, 1])
  
  # SE for rg using delta method: SE(rg) ≈ SE(h2_AB) / sqrt(h2_1 * h2_2)
  se_rg <- NA_real_
  if (!is.na(h2_1) && !is.na(h2_2) && h2_1 > 0 && h2_2 > 0) {
    se_rg <- se_h2_AB / sqrt(h2_1 * h2_2)
  }
  
  return(list(
    rg = rg,
    se_rg = se_rg,
    h2_AB = h2_AB,
    se_h2_AB = se_h2_AB,
    h2_1 = h2_1,
    h2_2 = h2_2,
    intercept = intercept,
    loglik = loglik_old,
    theta = theta,
    converged = converged,
    iterations = iter,
    mu = mu
  ))
}