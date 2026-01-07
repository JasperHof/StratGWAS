#' LD score regression
#'
#' Compute genetic covariance of subgroups using LD score regression
#'
#' @useDynLib StratGWAS, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#'
#' @param strata An object returned from stratify()
#' @param filename Prefix of genotype .bed file
#' @param nr_blocks Block size for reading in genotype data (default: 1000)
#' @param outfile Name of output file (should match previous step)
#' @param ss_list Optional: list of length K+2, containing input summary statistics specific to strata, binary trait, and stratification variable
#' @param lds Optional: data frame containing LD scores 
#' @return Returns covariance matrix of the strata
#' @export
gencov_ldsc <- function(strata, filename, nr_blocks = 1000, outfile, ss_list = NULL, lds = NULL) {

  # <performs some checks here> #

  # Read in data as multivariate phenotype
  K <- strata$K
  K_tot <- K + 2
  ids <- as.character(read.table(paste0(filename, ".fam"))[, 2])
  multi <- matrix(0, length(ids), strata$K)  
  for(k in 1:K) multi[, k] = strata[[paste0("group", k)]][match(ids,strata[[paste0("group", k)]][,1]), 3]

  # Scale phenotypes - but keep missings as missings
  for(k in 1:K){
    multi[!is.na(multi[,k]),k] = scale(as.numeric(multi[!is.na(multi[,k]),k]))
  }
  rownames(multi) <- ids

  # Perform a linear regression on the data
  #linear_gwas(filename, multi, nr_blocks, outfile)

  # Perform a linear regression on the data: (i) subtypes; (ii) disease; (iii) stratification variable
  if(is.null(ss_list)){
  multi = cbind(multi, as.numeric(scale(strata$y[match(ids, strata$y[,1]),3])), as.numeric(scale(strata$info[match(ids, strata$info[,1]),3])))
  multi_pheno <- cbind(ids, ids, multi)
  write.table(multi_pheno, paste0(outfile, ".strata"), quote = F, row = F, col = F)

  linear_gwas(filename, multi, nr_blocks, outfile)

  # Read in the linear regressions in list
  ss_list = vector("list", K_tot)
  for(k in 1:K_tot){
    ss_list[[k]] = read.table(paste0(outfile, ".pheno",k), head = T)
  }
  }

  # Compute LD scores for a subset of at most 1000 samples
  if(is.null(lds)){
  if (length(ids) > 1000) {
    geno_set <- sample(1:length(ids), size = 1000)
  } else {
    geno_set <- 1:length(ids)
  }
  geno_set = geno_set[order(geno_set)]

  lds = computeLDscoresFromBED(filename, geno_set)
  write.table(lds, paste0(outfile,".ldscores"), quote = F, row = F, col = T)
  }

  # Use LDSC for genetic correlations
  hers = rep(NA, K_tot)
  gencor = matrix(NA, K_tot, K_tot)
  lds_matched = lds$Tagging[match(ss_list[[1]]$Predictor, lds$Predictor)]

  for(i in 1:K_tot){
    for(j in i:K_tot){
      if(i == j){
        ldsc <- ldsc_cor(ss_list[[i]], ss_list[[j]], lds_matched)
        gencor[i,j] <- 1
        hers[i] <- max(as.numeric(ldsc$h2_1), 0.001)
      } else {
        ss1 = ss_list[[i]]
        ss2 = ss_list[[j]]
        ldsc = ldsc_cor(ss_list[[i]], ss_list[[j]], lds_matched)
        gencor[i,j] = gencor[j,i] = ldsc$rg
      }
    }
  }

  # Compute genetic covariance matrix
  gencov <- gencor
  for(k in 1:K_tot) gencov[k,] <- gencov[k,] * sqrt(hers[k])
  for(k in 1:K_tot) gencov[,k] <- gencov[,k] * sqrt(hers[k])

  write.table(hers, paste0(outfile,".hers"), quote = F, row = F, col = F)
  write.table(gencov, paste0(outfile,".gencov"), quote = F, row = F, col = F)
  write.table(gencor, paste0(outfile,".gencor"), quote = F, row = F, col = F)

  return(gencov)
}