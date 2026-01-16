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
stratgwas <- function(pheno, strat, filename, cov = NULL, block_size = 500) {
  
  # <performs some checks here> #

  fam <- read.table(paste0(filename, ".fam"))

  # reduce phenotype
  pheno <- pheno[pheno[, 1] %in% fam[, 1], ]
  ids <- pheno[,1]

  # reduce and scale stratification variable 
  strat <- strat[strat[, 1] %in% pheno[which(pheno[, 3] == 1), 1], ]
  #strat[,3] <- as.numeric(scale(rank(as.numeric(strat[, 3]))))
  strat[,3] <- as.numeric(scale(as.numeric(strat[, 3])))

  # create multivariate phenotype file for HE regression
  multi <- cbind(pheno[,3], 0, 0, 0)
  multi[match(strat[,1], ids),2] <- 1 + .5 * strat[, 3]
  multi[match(strat[,1], ids),3] <- 1 + .5 * strat[, 3]^2
  multi[match(strat[,1], ids),4] <- 1 + .5 * strat[, 3]^3
  rownames(multi) = ids

  # regress covariates from phenotype
  if(!is.null(cov)){
    cov_df <- cov[match(ids, cov[,1]), -(1:2)]
    cov_df[] <- lapply(cov_df, function(x) {
        x <- as.numeric(x)
        x[is.na(x)] <- mean(x, na.rm = TRUE)
        x
    })

    y <- pheno[match(ids, pheno[,1]), 3]
    names(y) <- ids
    rownames(cov_df) <- ids

    fit <- lm(y ~ ., data = cov_df)
    multi[match(names(residuals(fit)), ids), 1] <- residuals(fit)
  }

  # normalize all columns
  for(i in 1:dim(multi)[2]) multi[, i] = as.numeric(scale(multi[, i]))

  # make sure the phenotype is numerical matrix; IDs provided as rownames
  cor_g = he_multi_part(filename, multi, block_size)

  # compute genetic correlation
  cor_total = cor(multi, use = "complete.obs")
  cor_e = cor_total - cor_g

  # perform double decomposition
  P = eigen(cor_e)                                                    # decompose cor_e = U E UT
  S_inv = P$vectors %*% diag(1/sqrt(P$values)) %*% t(P$vectors)       # get S_inv = U E^{-1/2} U^T
  S = P$vectors %*% diag(sqrt(P$values)) %*% t(P$vectors)             # get S = U E^{1/2} U^T
  M = S_inv %*% as.matrix(cor_g) %*% S_inv                            # transform G into environment space - and decompose after
  V = eigen(M)$vectors                                                # now compute U for downstream use
  U1 = t(V) %*% S_inv

  # retrieve weightings and compute transformed phenotype
  weights = U1[1,]
  trans_pheno = cbind(ids, ids, multi %*% weights)

  # compute genetic correlation
  hers <- diag(cor_g)
  rg <- cor_g
  
  for(k in 1:dim(rg)[1]){
    if(hers[k] < 0){
      cat(sprintf("SNP heritability of trait %d is negative, so will not compute genetic correlation \n", k))
      rg[k,] <- NA
      rg[,k] <- NA
    } else {
      rg[k,] <- rg[k,] / sqrt(hers[k])
      rg[,k] <- rg[,k] / sqrt(hers[k])
    }
  }

  # compute inflation criterion
  gamma = rg[1,2]

  if(is.na(gamma)){
    message("Can not compute inflation criterion due to negative heritability.")
  } else {
    h2_Z = cor_g[2,2]

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