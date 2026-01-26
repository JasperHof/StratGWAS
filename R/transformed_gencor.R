#' Compute genetic correlation between transformed and original binary trait
#'
#' Given the genetic covariance matrix and transformation weights, compute
#' the genetic correlation between the transformed phenotype Y' = w^T Y
#' and the original binary trait Y[1].
#'
#' @param cor_g Genetic covariance matrix
#' @param weights Weight vector for transformation
#' @param binary_idx Index of binary trait (default = 1)
#'
#' @return List containing:
#'   - rg: Genetic correlation between transformed and original trait
#'   - h2_transformed: Heritability of transformed trait
#'   - h2_binary: Heritability of binary trait
#'   - cov_g: Genetic covariance between transformed and binary trait
#'
#' @details
#' For transformed phenotype Y' = w^T Y, the genetic correlation with Y[1] is:
#'   rg = Cov_g(Y', Y[1]) / sqrt(Var_g(Y') * Var_g(Y[1]))
#'
#' where:
#'   Cov_g(Y', Y[1]) = w^T * Sigma_g * e_1
#'   Var_g(Y') = w^T * Sigma_g * w
#'   Var_g(Y[1]) = Sigma_g[1,1]
#'
#' @export
compute_rg_transformed <- function(cor_g, weights, binary_idx = 1) {
  
  n_vars <- length(weights)
  
  # Create indicator vector for binary trait
  e1 <- numeric(n_vars)
  e1[binary_idx] <- 1
  
  # Genetic covariance between transformed and binary trait
  # Cov_g(Y', Y[1]) = w^T * Sigma_g * e_1
  cov_g <- as.numeric(t(weights) %*% cor_g %*% e1)
  
  # Variance of transformed trait
  # Var_g(Y') = w^T * Sigma_g * w
  var_g_transformed <- as.numeric(t(weights) %*% cor_g %*% weights)
  
  # Variance of binary trait
  var_g_binary <- cor_g[binary_idx, binary_idx]
  
  # Genetic correlation
  rg <- cov_g / sqrt(var_g_transformed * var_g_binary)
  
  return(list(
    rg = rg,
    h2_transformed = var_g_transformed,
    h2_binary = var_g_binary,
    cov_g = cov_g
  ))
}


#' Compute genetic correlation matrix for multiple transformed variables
#'
#' @param cor_g Genetic covariance matrix
#' @param weights_matrix Matrix where each column is a weight vector
#'
#' @return Genetic correlation matrix including original and transformed variables
#'
#' @export
compute_rg_matrix <- function(cor_g, weights_matrix) {
  
  n_orig <- nrow(cor_g)
  n_transformed <- ncol(weights_matrix)
  n_total <- n_orig + n_transformed
  
  # Initialize correlation matrix
  rg_matrix <- matrix(NA, n_total, n_total)
  
  # Fill in original variable correlations
  # First convert covariances to correlations
  h_orig <- sqrt(diag(cor_g))
  rg_orig <- cor_g / outer(h_orig, h_orig)
  rg_matrix[1:n_orig, 1:n_orig] <- rg_orig
  
  # Compute heritabilities of transformed variables
  h_transformed <- numeric(n_transformed)
  for(i in 1:n_transformed) {
    w <- weights_matrix[, i]
    h_transformed[i] <- sqrt(as.numeric(t(w) %*% cor_g %*% w))
  }
  
  # Compute correlations between original and transformed
  for(i in 1:n_orig) {
    e_i <- numeric(n_orig)
    e_i[i] <- 1
    
    for(j in 1:n_transformed) {
      w_j <- weights_matrix[, j]
      cov_ij <- as.numeric(t(w_j) %*% cor_g %*% e_i)
      rg_ij <- cov_ij / (h_orig[i] * h_transformed[j])
      
      rg_matrix[i, n_orig + j] <- rg_ij
      rg_matrix[n_orig + j, i] <- rg_ij
    }
  }
  
  # Compute correlations between transformed variables
  for(i in 1:n_transformed) {
    for(j in 1:n_transformed) {
      w_i <- weights_matrix[, i]
      w_j <- weights_matrix[, j]
      
      cov_ij <- as.numeric(t(w_i) %*% cor_g %*% w_j)
      rg_ij <- cov_ij / (h_transformed[i] * h_transformed[j])
      
      rg_matrix[n_orig + i, n_orig + j] <- rg_ij
    }
  }
  
  return(rg_matrix)
}


#' Compute expected improvement in GWAS power
#'
#' Estimate the expected improvement in GWAS power from using transformed trait
#' based on the genetic correlation and heritability changes.
#'
#' @param h2_original Heritability of original trait
#' @param h2_transformed Heritability of transformed trait
#' @param rg Genetic correlation between original and transformed
#' @param n_samples Sample size
#'
#' @return List with power improvement metrics
#'
#' @details
#' The effective sample size for detecting a variant is approximately:
#'   N_eff = N * h2 * rg^2
#'
#' So the improvement factor is:
#'   improvement = (h2_transformed * rg^2) / h2_original
#'
#' @export
compute_power_improvement <- function(h2_original, h2_transformed, rg, n_samples = NULL) {
  
  # Improvement in effective sample size
  improvement_factor <- (h2_transformed * rg^2) / h2_original
  
  # Approximate improvement in chi-square statistic
  # chi2 ~ N * h2 * effect^2
  chi2_improvement <- improvement_factor
  
  # Approximate improvement in z-score (sqrt of chi2)
  z_improvement <- sqrt(improvement_factor)
  
  result <- list(
    improvement_factor = improvement_factor,
    chi2_ratio = chi2_improvement,
    z_ratio = z_improvement,
    expected_power_gain = (improvement_factor - 1) * 100  # percentage
  )
  
  if(!is.null(n_samples)) {
    result$n_eff_original <- n_samples * h2_original
    result$n_eff_transformed <- n_samples * h2_transformed * rg^2
  }
  
  return(result)
}


#' Diagnostic plot for transformed trait
#'
#' @param cor_g Genetic covariance matrix
#' @param weights Weight vector
#' @param var_names Optional variable names
#'
#' @export
plot_transformation_diagnostics <- function(cor_g, weights, var_names = NULL) {
  
  if(is.null(var_names)) {
    var_names <- paste0("Var", 1:length(weights))
  }
  
  # Compute genetic correlation
  rg_result <- compute_rg_transformed(cor_g, weights)
  
  # Create summary
  cat("=== Transformation Diagnostics ===\n\n")
  cat(sprintf("Original trait heritability:    %.4f\n", rg_result$h2_binary))
  cat(sprintf("Transformed trait heritability: %.4f\n", rg_result$h2_transformed))
  cat(sprintf("Genetic correlation (rg):       %.4f\n\n", rg_result$rg))
  
  # Power improvement
  power <- compute_power_improvement(
    rg_result$h2_binary,
    rg_result$h2_transformed,
    rg_result$rg
  )
  
  cat("=== Expected GWAS Power Improvement ===\n\n")
  cat(sprintf("Improvement factor:  %.2fx\n", power$improvement_factor))
  cat(sprintf("Power gain:          %.1f%%\n\n", power$expected_power_gain))
  
  # Weight contributions
  cat("=== Variable Weights ===\n\n")
  weight_df <- data.frame(
    Variable = var_names,
    Weight = round(weights, 4),
    AbsWeight = round(abs(weights), 4)
  )
  weight_df <- weight_df[order(-weight_df$AbsWeight), ]
  print(weight_df, row.names = FALSE)
  
  cat("\n")
  
  # Warning if rg is too low
  if(abs(rg_result$rg) < 0.7) {
    warning("Genetic correlation < 0.7 suggests transformed trait may not be ",
            "capturing the same genetic signal as the original trait!")
  }
  
  invisible(list(
    rg_result = rg_result,
    power = power,
    weights = weight_df
  ))
}