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

  # Convert to data.frame
  pheno <- as.data.frame(pheno)
  strat <- as.data.frame(strat)

  # Store IDs
  ids <- pheno[, 1]

  # Identify cases with missing stratification variables
  cases <- pheno[which(pheno[, 3] == 1), 1]
  cases_nostrat <- strat[which(strat[, 1] %in% cases & is.na(strat[, 3])), 1]

  # These are removed for now and later added
  pheno <- pheno[!(pheno[, 1] %in% cases_nostrat), ]
  strat <- strat[!(strat[, 1] %in% cases_nostrat), ]

  # Match with phenotype
  strat <- strat[match(pheno[, 1], strat[, 1]), ]

  # Extract stratification variable and compute quintiles
  strat_cases <- strat[which(pheno[,3] == 1), ]

  # Determine stratification approach
  n_unique <- length(unique(strat[!is.na(strat[, 3]), 3]))
  sparse <- n_unique < 10

  # Define groups 
  if (sparse) {
    # Groups defined by input groups
    message(sprintf("Only %d unique stratification values detected. 
                    Performing sparse StratGWAS with K=%d strata",
                    n_unique, n_unique))
    K <- n_unique
    strat_cases$groups <- match(strat_cases[, 3], sort(unique(strat[, 3])))
  } else {
    # Quantile-based stratification with jittered ranks for tie-breaking
    strat_cases$groups <- assign_to_quantiles(strat_cases[, 3], K)
  }

  # Create stratified phenotype lists
  strata <- create_strata_list(pheno, strat_cases, K)

  # Return list with information
  strata[["K"]] = K
  strata[["y"]] = pheno
  strata[["Z"]] = strat
  strata[["info"]] = strat_cases
  strata[["ids"]] = ids
  strata[["strat_miss"]] = cases_nostrat
  strata[["sparse"]] = sparse

  return(strata)
}

#' Assign values to quantile-based groups
#'
#' @param x Numeric vector to stratify
#' @param K Number of quantile groups
#' @return Integer vector of group assignments (1 to K)
assign_to_quantiles <- function(x, K = 5) {

  # Calculate ranks (handles ties by averaging)
  ranks <- rank(x, ties.method = "average")

  # Convert to percentiles
  percentiles <- ranks / length(ranks)

  # Add minimal jitter to break remaining ties
  set.seed(123)  # For reproducibility
  percentiles_jittered <- percentiles + rnorm(length(percentiles),
                                              mean = 0,
                                              sd = 1e-6)

  # Assign to quantile groups
  breaks <- seq(0, 1, length.out = K + 1)
  groups <- cut(percentiles_jittered,
                breaks = breaks,
                labels = 1:K,
                include.lowest = TRUE)

  return(as.integer(groups))
}

#' Create list of stratified phenotypes
#'
#' @param pheno Full phenotype data.frame
#' @param strat_cases Stratification info for cases
#' @param K Number of strata
#' @return List of phenotype data.frames, one per stratum
create_strata_list <- function(pheno, strat_cases, K) {

  strata <- vector("list", K)
  names(strata) <- paste0("group", 1:K)

  for (k in 1:K) {
    # Select controls (pheno == 0) and cases from stratum k
    cases_in_stratum <- strat_cases[strat_cases$groups == k, 1]
    stratum_pheno <- pheno[pheno[, 3] == 0 | pheno[, 1] %in% cases_in_stratum, ]

    # Normalize phenotype (mean 0, sd 1)
    stratum_pheno[, 3] <- as.numeric(scale(stratum_pheno[, 3]))

    strata[[k]] <- stratum_pheno
  }

  return(strata)
}