#' Compute transformed phenotype using StratGWAS
#'
#' Compute genetic covariance of subgroups using Haseman-Elston regression
#'
#' @useDynLib StratGWAS, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' 
#' @param pheno Matrix containing binary outcome in PLINK format
#' @param strat Matrix containing stratification variable in PLINK format
#' @param filename Prefix of input .bed file
#' @return Returns covariance matrix of the strata
#' @export
stratgwas <- function(pheno, strat, filename, cov = NULL, block_size = 500, cor_g = NULL, strat_mat = NULL, alpha = 0) {
  
  # <performs some checks here> #

  fam <- read.table(paste0(filename, ".fam"))

  # reduce phenotype
  pheno <- pheno[pheno[, 1] %in% fam[, 1], ]
  ids <- pheno[,1]

  if(is.null(strat_mat)){
    
    # reducing stratification variable
    strat <- strat[strat[, 1] %in% pheno[which(pheno[, 3] == 1), 1], ]
    strat[is.na(strat[, 3]), 3] <- mean(strat[is.na(strat[, 3]), 3], na.rm = TRUE)
    strat[,3] <- as.numeric(scale(as.numeric(strat[, 3])))

    # construct multivariate phenotype
    strat_mat <- spline_manual(strat[, 3])

    #strat_cov <- strat[match(ids, strat[, 1]), 3]
    #strat_cov[is.na(strat_cov)] <- mean(strat_cov, na.rm = T)
    #multi <- cbind(pheno[, 3], strat_cov, matrix(0, nrow = nrow(pheno), ncol = ncol(strat_mat)))

    multi <- cbind(pheno[, 3], matrix(0, nrow = nrow(pheno), ncol = ncol(strat_mat)))
    multi[match(strat[, 1], ids), -1] <- alpha + strat_mat

    rownames(multi) <- ids
  } else {
    multi <- strat_mat
    rownames(multi) <- ids
  }

  # create multivariate phenotype file for HE regression
  #multi <- cbind(pheno[,3], 0, 0, 0, 0)
  #multi[match(strat[,1], ids),2] <- 0.5 + .5 * strat[, 3]
  #multi[match(strat[,1], ids),3] <- 0.5 + .5 * strat[, 3]^2
  #multi[match(strat[,1], ids),4] <- 0.5 + .5 * strat[, 3]^3
  #multi[match(strat[,1], ids),4] <- 0.5 + .5 * sin(strat[, 3])
  #multi[match(strat[,1], ids),5] <- 0.5 + .5 * cos(strat[, 3])
  #rownames(multi) <- ids

  # regress covariates from phenotype
  if(!is.null(cov)){
    cov_df <- cov[match(ids, cov[,1]), -(1:2), drop = F]
    cov_df[] <- lapply(cov_df, function(x) {
        x <- as.numeric(x)
        x[is.na(x)] <- mean(x, na.rm = TRUE)
        x
    })
    rownames(cov_df) <- ids

    # regress from all columns of multi
    for(i in 1:ncol(multi)){
      y <- multi[, i]
      names(y) <- ids

      fit <- lm(y ~ ., data = cov_df)
      multi[match(names(residuals(fit)), ids), i] <- residuals(fit)
    }
  }

  # normalize all columns
  for(i in 1:dim(multi)[2]) multi[, i] = as.numeric(scale(multi[, i]))

  # compute genetic covariance using randomized HE regression
  if(is.null(cor_g)) cor_g = he_multi_part(filename, multi, block_size)

  ############
  # compute genetic correlation
  ############

  rg <- cor_g

  # check for negative heritability
  hers <- diag(cor_g)
  her_neg <- which(hers < 0)
  
  for(k in 1:dim(rg)[1]){
    if(hers[k] < 0) {
      cat(sprintf("SNP heritability of trait %d is negative, so will not compute genetic correlation \n", k))
      rg[k, ] <- NA
      rg[, k] <- NA
    } else {
      rg[k, ] <- rg[k, ] / sqrt(hers[k])
      rg[, k] <- rg[, k] / sqrt(hers[k])
    }
  }

  # store genetic correlation between binary trait and auxiliary variable
  gamma <- rg[1, 2]
  h2_Z <- cor_g[2, 2]

  # now remove this variable from cor_g and multi, not further used
  #cor_g <- cor_g[-2, -2]
  #multi <- multi[, -2]

  # compute environmenta; correlation
  cor_total = cor(multi, use = "complete.obs")
  cor_e = cor_total - cor_g

  if(length(her_neg) > 0){
    # Create informative message about which terms are being removed
    term_names <- c("binary", "linear", "quadratic", "cubic")[her_neg]
    message("Heritability of ", paste(term_names, collapse = ", "),
            " term(s) is negative and will be ignored.")

    if(length(her_neg) >= 3){
      warning("Too many terms have negative heritability.")
    }

    # Remove negative heritability terms
    cor_g_use <- cor_g[-her_neg, -her_neg]
    cor_e_use <- cor_e[-her_neg, -her_neg]
    multi_use <- multi[, -her_neg]
  } else {
    # Use full matrices if no negative heritabilities
    cor_g_use <- cor_g
    cor_e_use <- cor_e
    multi_use <- multi
  }

  # perform double decomposition
  P <- eigen(cor_e_use)
  S_inv <- P$vectors %*% diag(1/sqrt(P$values)) %*% t(P$vectors)
  S <- P$vectors %*% diag(sqrt(P$values)) %*% t(P$vectors)
  M <- S_inv %*% as.matrix(cor_g_use) %*% S_inv
  V <- eigen(M)$vectors
  U1 <- t(V) %*% S_inv

  # retrieve weightings and compute transformed phenotype
  weights <- U1[1,]
  trans_pheno <- cbind(ids, ids, multi_use %*% weights)

  # compute inflation criterion
  if(is.na(gamma)){
    message("Can not compute inflation criterion due to negative heritability.")
  } else {
    # compute a2 - the correlation between stratification Y' and Z
    y_trans = as.numeric(scale(trans_pheno[,3]))
    y = as.numeric(scale(multi[,1]))
    z = as.numeric(scale(multi[,2]))

    # fit regression
    fit = lm(y_trans ~ y + z - 1)
    coefs = fit$coefficients
    a2 = coefs["z"]

    # Compute inflation criterion using genetic correlation
    criterion = a2^2 * (1 - gamma^2) * h2_Z
    message(paste0("Expected inflation criterion is ",round(criterion, 4)))
    if(criterion > 0.01) warning("Inflation criterion is greater than 0.01.")
  }

  # Return list with information
  object = vector("list")
  object[["pheno"]] = trans_pheno
  object[["cor_g"]] = cor_g
  object[["cor_e"]] = cor_e
  object[["rg"]] = rg
  object[["criterion"]] = criterion
  object[["weights"]] = weights

  return(object)
}