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
  strat[,3] <- as.numeric(scale(rank(as.numeric(strat[, 3]))))

  # create multivariate phenotype file for HE regression
  multi <- cbind(pheno[,3], 0, 0, 0)
  multi[match(strat[,1], ids),2] <- 1 + .5 * strat[, 3]
  multi[match(strat[,1], ids),3] <- 1 + .5 * strat[, 3]^2
  multi[match(strat[,1], ids),4] <- 1 + .5 * strat[, 3]^3
  rownames(multi) = ids

  # regress covariates from phenotype
  if(!is.null(cov)){
    
    # create covariate matrix for performing linear regression
    cov_mat <- cbind(ids, ids, cov[match(ids, cov[, 1]), -c(1, 2)])
    for(i in 3:dim(cov_mat)[2]) cov_mat[,i] <- as.numeric(cov_mat[,i])
    for(i in 3:dim(cov_mat)[2]) cov_mat[,i] <- mean(cov_mat[,i], na.rm = T)
    
    y <- pheno[, 3]
    names(y) <- ids

    fit <- lm(y ~ cov_mat)
    multi[match(names(fit$residuals), ids), 1] <- fit$residuals
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
  #U1 %*% t(cor_e) %*% t(U1)
  #U1 %*% t(cor_g) %*% t(U1)

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

  # Return list with information
  object = vector("list")
  object[["pheno"]] = trans_pheno
  object[["cor_g"]] = cor_g
  object[["cor_e"]] = cor_e
  object[["rg"]] = rg
  object[["weights"]] = weights

  return(object)
}