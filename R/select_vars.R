#' Greedy forward selection of variables for heritability improvement
#'
#' @param cor_g Genetic covariance matrix (full)
#' @param cor_e Environmental covariance matrix (full)
#' @param min_improvement Minimum relative improvement required (default 0.05 = 5%)
#' @param verbose Print progress information
#'
#' @return List containing:
#'   - selected_vars: Indices of selected variables (including binary trait at position 1)
#'   - h2_trajectory: Vector of heritabilities at each step
#'   - improvements: Vector of relative improvements at each step
#'   - weights: Final weight vector for selected variables
#'
#' @keywords internal
select_vars <- function(cor_g, cor_e, min_improvement = 0.05, verbose = TRUE) {
  
  n_vars <- nrow(cor_g)
  
  # Start with just the binary trait (index 1)
  selected <- c(1)
  candidates <- setdiff(2:n_vars, selected)
  
  # Compute baseline heritability with just binary trait
  h2_current <- compute_transformed_h2(cor_g, cor_e, selected)
  
  h2_trajectory <- h2_current
  improvements <- numeric(0)
  
  if(verbose) {
    cat(sprintf("Starting heritability (binary trait only): %.4f\n\n", h2_current))
  }
  
  # Greedy forward selection
  while(length(candidates) > 0) {
    
    best_h2 <- h2_current
    best_var <- NULL
    best_improvement <- 0
    
    # Try adding each candidate variable
    for(var in candidates) {
      test_vars <- c(selected, var)
      h2_test <- compute_transformed_h2(cor_g, cor_e, test_vars)
      
      # Compute relative improvement
      improvement <- (h2_test - h2_current) / h2_current
      
      if(h2_test > best_h2) {
        best_h2 <- h2_test
        best_var <- var
        best_improvement <- improvement
      }
    }
    
    # Check if best improvement meets threshold
    if(best_improvement >= min_improvement) {
      selected <- c(selected, best_var)
      candidates <- setdiff(candidates, best_var)
      h2_current <- best_h2
      
      h2_trajectory <- c(h2_trajectory, h2_current)
      improvements <- c(improvements, best_improvement)
      
      if(verbose) {
        cat(sprintf("Added variable %d: h2 = %.4f (improvement: %.2f%%)\n", 
                    best_var, h2_current, best_improvement * 100))
      }
    } else {
      # No variable meets improvement threshold, stop
      if(verbose) {
        cat(sprintf("\nStopping: best improvement (%.2f%%) below threshold (%.2f%%)\n",
                    best_improvement * 100, min_improvement * 100))
      }
      break
    }
  }
  
  if(verbose) {
    cat(sprintf("\nFinal heritability: %.4f\n", h2_current))
    cat(sprintf("Selected %d/%d variables\n", length(selected), n_vars))
  }
  
  # Compute final weights for selected variables
  cor_g_selected <- cor_g[selected, selected, drop = FALSE]
  cor_e_selected <- cor_e[selected, selected, drop = FALSE]
  
  weights_result <- compute_weights(cor_g_selected, cor_e_selected)
  
  # Create full weight vector (zeros for unselected variables)
  weights_full <- numeric(n_vars)
  weights_full[selected] <- weights_result$weights
  
  return(list(
    selected_vars = selected,
    h2_trajectory = h2_trajectory,
    improvements = improvements,
    weights = weights_full,
    weights_selected = weights_result$weights,
    h2_final = h2_current
  ))
}


#' Compute transformed heritability for a subset of variables
#'
#' @param cor_g Genetic covariance matrix (full)
#' @param cor_e Environmental covariance matrix (full)
#' @param vars Indices of variables to include
#'
#' @return Heritability of the first principal component
#'
#' @keywords internal
compute_transformed_h2 <- function(cor_g, cor_e, vars) {
  
  # Subset matrices
  cor_g_sub <- cor_g[vars, vars, drop = FALSE]
  cor_e_sub <- cor_e[vars, vars, drop = FALSE]
  
  # Handle case with single variable
  if(length(vars) == 1) {
    return(cor_g_sub[1, 1])
  }
  
  # Eigendecomposition of environmental covariance
  P <- eigen(cor_e_sub)
  
  # Check for numerical issues
  if(any(P$values <= 0)) {
    warning("Non-positive eigenvalues in environmental covariance")
    P$values[P$values <= 0] <- 1e-8
  }
  
  S_inv <- P$vectors %*% diag(1/sqrt(P$values)) %*% t(P$vectors)
  
  # Transform genetic covariance
  M <- S_inv %*% as.matrix(cor_g_sub) %*% S_inv
  
  # Eigendecomposition to find optimal linear combination
  V <- eigen(M)$vectors
  U1 <- t(V) %*% S_inv
  
  # Compute heritability of transformed phenotype
  D <- U1 %*% cor_g_sub %*% t(U1)
  
  return(as.numeric(D[1, 1]))
}


#' Compute weights for selected variables
#'
#' @param cor_g_sub Genetic covariance matrix (subset)
#' @param cor_e_sub Environmental covariance matrix (subset)
#'
#' @return List with weights and other decomposition results
#'
#' @keywords internal
compute_weights <- function(cor_g_sub, cor_e_sub) {
  
  # Eigendecomposition of environmental covariance
  P <- eigen(cor_e_sub)
  
  if(any(P$values <= 0)) {
    P$values[P$values <= 0] <- 1e-8
  }
  
  S_inv <- P$vectors %*% diag(1/sqrt(P$values)) %*% t(P$vectors)
  
  # Transform genetic covariance
  M <- S_inv %*% as.matrix(cor_g_sub) %*% S_inv
  
  # Eigendecomposition
  V <- eigen(M)$vectors
  U1 <- t(V) %*% S_inv
  
  # First row of U1 gives the weights
  weights <- U1[1, ]
  
  # Compute heritability
  D <- U1 %*% cor_g_sub %*% t(U1)
  h2 <- as.numeric(D[1, 1])
  
  return(list(
    weights = weights,
    h2 = h2,
    U1 = U1,
    M = M
  ))
}