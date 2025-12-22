#' LD score regression
#'
#' Compute genetic covariance of subgroups using LD score regression
#'
#' @useDynLib StratGWAS, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' 
#' @param strata An object returned from stratify()
#' @param filename Prefix of genotype .bed file
#' @return Returns covariance matrix of the strata
#' @export
gencov_ldsc <- function(strata, filename, nr_blocks = 1000, outfile, tag) {
  
  # <performs some checks here> #

  bim <- read.table(paste0(filename, ".bim"))
  overlap <- intersect(bim[,2], tag$Predictor)

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
  cor = linear_gwas(filename, multi, nr_blocks, outfile)

  # Read in the linear regressions in list
  ss_list = vector("list", K)
  for(k in 1:K){
    ss = read.table(paste0(outfile, ".pheno",k), head = T)
    ss_list[[k]] = ss[match(overlap, ss$Predictor),]
  }

  gencor = matrix(NA, K, K)
  for(i in 1:K){
    for(j in i:K){
      if(i == j){
        gencor[i,j] = ldsc(ss_list[[i]], tag$Tagging[match(overlap, tag$Predictor)])
      } else {
        cor = ldsc(ss_list[[i]], ss_list[[j]], tag$Tagging[match(overlap, tag$Predictor)])
        gencor[i,j] = gencor[j,i] = cor$rg * sqrt(cor$h2_1 * cor$h2_2)
      }
    }
  }

  return(gencor)
}