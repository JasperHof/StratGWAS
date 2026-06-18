#' Compute genetic liability for all individuals
#'
#' !!! Note that this function is in development
#' 
#' @export
liability <- function(pheno, filename, strat_cont = NULL, strat_cat = NULL, env = T) {

  # Create multivariate phenotype
  multi <- pheno
  colnames(multi) <- c("FID", "IID", "Pheno")

  if (!is.null(strat_cont)) {
    for (j in 3:ncol(strat_cont)){
      cont_mean <- mean(strat_cont[match(pheno[,2], strat_cont[, 2]), j], na.rm = T)

      multi <- cbind(multi, strat_cont[match(pheno[,2], strat_cont[, 2]), j], 
              (strat_cont[match(pheno[,2], strat_cont[, 2]), j] - cont_mean)^2)

      colnames(multi)[(ncol(multi) - 1):ncol(multi)] <- paste0("cont_", j-2, "_", 1:2)
    }
  }
  
  if (!is.null(strat_cat)) {
    for (j in 3:ncol(strat_cat)){
      multi <- cbind(multi, strat_cat[match(pheno[,2], strat_cont[, 2]), j])
      colnames(multi)[ncol(multi)] <- paste0("cat_", j-2)
    }
  }

  # Impute means and normalize
  for (j in 3:ncol(multi)) {
    multi[is.na(multi[, j]), j] <- mean(multi[, j], na.rm = T)
    multi[, j] <- as.numeric(scale(multi[, j]))
  }

  # HE regression on the phenotype
  multi_he <- as.matrix(multi[, -c(1,2)])
  rownames(multi_he) <- multi[, 2]

  gencov <- he_multi_part(filename, multi_he, 1000)

  # Should I scale to liability? Don't think so

  # Compute adjusted covariance matrix
  gencov_adj <- gencov

  for (i in 2:nrow(gencov_adj)) {
    gencov_adj[i, -1] <- gencov_adj[i, -1] * sqrt(abs(gencov_adj[i, 1] / gencov_adj[i, i]))
    gencov_adj[-1, i] <- gencov_adj[-1, i] * sqrt(abs(gencov_adj[i, 1] / gencov_adj[i, i]))
  }

  # Adjust to positive definite matrix
  gencov_adj <- as.matrix(nearPD(gencov_adj)$mat)

  # Get total and environmental correlations
  if (env == T) {
    total_cor <- cor(multi_he, use = 'complete.obs')
    env_cor <- total_cor - gencov

    # Adjust to positive definite matrix
    env_cor <- as.matrix(nearPD(env_cor)$mat)

    # Do a double eigendecomposition
    Wmat = eigen(env_cor)                                                  # decompose E = U diag(Le) UT
    Wmat4 = Wmat$vectors %*% diag(1/sqrt(Wmat$values)) %*% t(Wmat$vectors) # get Wmat4 = U E^{-1/2} U^T
    Wmat5 = Wmat$vectors %*% diag(sqrt(Wmat$values)) %*% t(Wmat$vectors)   # get Wmat5 = U E^{1/2} U^T
    M = Wmat3 = Wmat4 %*% as.matrix(gencov_adj) %*% Wmat4               # transform G into environment space (???) - and decompose after
    V = eigen(Wmat3)$vectors                                               # now compute U for downstream use
    U1 = t(V) %*% Wmat4
    U2 = Wmat5 %*% V

    # does it work both for covs and cors_env?
    #U1 %*% t(env_cor) %*% t(U1)
    #U1 %*% t(gencov_adj) %*% t(U1)

    weights <- U1[1, ]
  } else {
    Wmat <- eigen(gencov_adj)
    weights <- Wmat$vectors[, 1]
  }
  

  # Compute new phenotype
  trans <- cbind(multi[, c(1, 2)], as.matrix(multi[, -c(1, 2)]) %*% weights)

  return(trans)
}