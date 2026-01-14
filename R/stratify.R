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

  # Track individuals with missing stratification variable
  cases <- pheno[which(pheno[, 3] == 1),1]
  cases_nostrat <- strat[which(strat[, 1] %in% cases & is.na(strat[, 3])), 1]

  # These are removed for now and later added
  pheno <- pheno[!pheno[,1] %in% cases_nostrat,]
  strat <- strat[!strat[,1] %in% cases_nostrat,]

  # Read in summary statistics
  strata <- vector("list")
  strat = as.data.frame(strat[match(pheno[,1], strat[,1]),]) # match stratification variable with phenotype

  # Extract stratification variable and compute quintiles
  strat_cases = strat[which(pheno[,3] == 1),]
  strat_cases$order = rank(strat_cases[,3])/length(strat_cases[,3])

  # Define groups 
  if(length(unique(strat[,3])) >= 10){
    strat_cases$groups = cut(strat_cases$order, breaks = c(0,1:(K-1)/K,Inf), labels = 1:5)
    sparse = F
  } else {
    strat_cases$groups = match(strat_cases[,3], unique(strat[,3]))
    K = length(unique(strat[!is.na(strat[,3]),3]))
    sparse = T
    message("Less than 10 unique stratification values, so will performe sparse StratGWAS")
  }

  # Create list with phenotype file for separate NORMALIZED strata
  for(k in 1:K){
    strata[[paste0("group",k)]] = pheno[which(pheno[,3] == 0 | pheno[,1] %in% strat_cases[strat_cases$groups == k,1]),]
    strata[[paste0("group",k)]][,3] = as.numeric(scale(strata[[paste0("group",k)]][,3]))
  }

  # Return list with information
  strata[["K"]] = K
  strata[["y"]] = pheno
  strata[["Z"]] = strat
  strata[["info"]] = strat_cases
  strata[["strat_miss"]] = cases_nostrat
  strata[["sparse"]] = sparse

  return(strata)
}