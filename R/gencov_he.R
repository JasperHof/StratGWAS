#' Haseman-Elston regression
#'
#' Compute genetic covariance of subgroups using Haseman-Elston regression
#'
#' @param strata An object returned from stratify()
#' @param filename Prefix of genotype .bed file
#' @return Returns covariance matrix of the strata
#' @export
gencov_he <- function(strata, filename) {
  
  # <performs some checks here> #

  bim = read.table(paste0(filename,".bim"))

  # Read in data as multivariate phenotype
  K = strata$K
  ids = strata$y[,1]
  multi = matrix(0, length(ids), strata$K)  
  for(k in 1:K) multi[,k] = strata[[paste0("group",k)]][match(ids,strata[[paste0("group",k)]][,1]),3]

  # Impute missing values and keep track of missingness
  miss = rep(0,K)
  for(k in 1:K){
    miss[k] = mean(is.na(multi[,k]))
    multi[is.na(multi[,k]),k] = mean(multi[,k], na.rm=T)
    multi[,k] = scale(as.numeric(multi[,k]))
  }
  rownames(multi) = ids

  # make sure the phenotype is numerical matrix; IDs provided as rownames
  cor = he_multi_part(filename, multi, 1000)

  return(cor)
}