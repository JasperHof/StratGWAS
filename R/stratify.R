#' Stratify cases
#'
#' Compute case subgroups by stratification based on input variable
#'
#' @param phenofile Binary input phenotype in PLINK format
#' @param stratfile Stratification variable in PLINK format
#' @param filename Prefix of input .bed file
#' @return Returns list containing subgroup phenotypes
#' @export
stratify <- function(pheno, strat, K = 5) {
  
  # Check input data
  stratify_checks(pheno, strat, K)

  # Read in summary statistics
  strata <- vector("list")
  strat = as.data.frame(strat[match(pheno[,1], strat[,1]),]) # match stratification variable with phenotype

  # <performs some checks here> #
  # what to do with missing stratification variable?

  # Extract stratification variable and compute quintiles
  strat_cases = strat[which(pheno[,3] == 1),]
  strat_cases[,3][is.na(strat_cases[,3])] = mean(strat_cases[,3], na.rm = T)
  strat_cases$order = rank(strat_cases[,3])/length(strat_cases[,3])
  strat_cases$groups = cut(strat_cases$order, breaks = c(0,1:(K-1)/K,Inf), labels = 1:5)

  # Create list with phenotype file for separate NORMALIZED strata
  for(k in 1:K){
    strata[[paste0("group",k)]] = pheno[which(pheno[,3] == 0 | pheno[,1] %in% strat_cases[strat_cases$groups == k,1]),]
    strata[[paste0("group",k)]][,3] = as.numeric(scale(strata[[paste0("group",k)]][,3]))
  }

  # < build in checks - check case distribution >
  strata[["K"]] = K
  strata[["y"]] = pheno
  strata[["Z"]] = strat
  strata[["info"]] = strat_cases

  return(strata)
}