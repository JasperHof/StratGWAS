#' LD score regression
#'
#' Compute genetic covariance of subgroups using LD score regression
#'
#' @useDynLib StratGWAS, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' 
#' @param strata An object returned from stratify()
#' @param filename Prefix of genotype .bed file
#' @param outfile Name of output file for strata-specific summary statistics
#' @param tag Numerical matrix containing LD scores
#' @return Returns covariance matrix of the strata
#' @export
gencov_ldsc <- function(strata, filename, nr_blocks = 1000, outfile) {
  
  # <performs some checks here> #

  # Read in data as multivariate phenotype
  K = strata$K
  ids = strata$y[,1]
  multi = matrix(0, length(ids), strata$K)  
  for(k in 1:K) multi[,k] = strata[[paste0("group",k)]][match(ids,strata[[paste0("group",k)]][,1]),3]

  # Scale phenotypes - but keep missings as missings
  for(k in 1:K){
    multi[!is.na(multi[,k]),k] = scale(as.numeric(multi[!is.na(multi[,k]),k]))
  }
  rownames(multi) = ids

  # Perform a linear regression on the data
  #linear_gwas(filename, multi, nr_blocks, outfile)

  # Perform a linear regression on the data: (i) subtypes; (ii) disease; (iii) stratification variable
  multi = cbind(multi, as.numeric(scale(strata$y[,3])), as.numeric(scale(strata$info[match(strata$y[,1], strata$info[,1]),3])))
  linear_gwas(filename, multi, nr_blocks, outfile)

  # Read in the linear regressions in list
  K_tot = K + 2
  ss_list = vector("list", K_tot)
  for(k in 1:K_tot){
    ss_list[[k]] = read.table(paste0(outfile, ".pheno",k), head = T)
  }

  # Compute LD scores for a subset of 1000 samples
  if (length(ids) > 1000) {
    geno_set <- sample(1:length(ids), size = 1000)
  } else {
    geno_set <- 1:length(ids)
  }
  geno_set = geno_set[order(geno_set)]

  lds = computeLDscoresFromBED(filename, geno_set)

  # Use for genetic correlations
  gencor = matrix(NA, K_tot, K_tot)
  lds_matched = lds$Tagging[match(ss_list[[1]]$Predictor, lds$Predictor)]

  for(i in 1:K_tot){
    for(j in i:K_tot){
      if(i == j){
        gencor[i,j] = ldsc(ss_list[[i]], lds_matched)
      } else {
        ss1 = ss_list[[i]]
        ss2 = ss_list[[j]]
        cor = ldsc_cor(ss_list[[i]], ss_list[[j]], lds_matched)
        gencor[i,j] = gencor[j,i] = cor$cov_g
      }
    }
  }

  return(gencor)
}