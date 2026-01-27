#' Stratify cases
#'
#' Compute case subgroups by stratification based on input variable
#'
#' @param phenofile Binary input phenotype in PLINK format
#' @param stratfile Stratification variable in PLINK format
#' @param filename Prefix of input .bed file
#' @return Returns list containing subgroup phenotypes
#' @export
stratify <- function(pheno, strat, K = 5, cov = NULL) {

  # Check input data
  stratify_checks(pheno, strat, K)

  # Convert to data.frame
  pheno <- as.data.frame(pheno)
  strat <- as.data.frame(strat)

  # Store IDs
  ids <- pheno[, 1]

  # Match with phenotype
  strat <- strat[match(ids, strat[, 1]), ]
  strat[, 1] <- strat[, 2] <- ids

  # Identify cases with missing stratification variables
  cases <- pheno[which(pheno[, 3] == 1), 1]
  cases_nostrat <- strat[which(strat[, 1] %in% cases & is.na(strat[, 3])), 1]

  # Extract stratification variable and compute quintiles
  strat_cases <- strat[which(pheno[, 3] == 1 & !(pheno[, 1] %in% cases_nostrat)), ]

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
    result <- assign_to_quantiles(strat_cases[, 3], K)
    strat_cases$groups <- result$groups
    strat_cases$order <- result$order
  }

  # Create stratified phenotype lists
  strata <- create_strata_list(pheno, strat_cases, K, cov = cov, ids = ids)

  # Return list with information
  strata[["K"]] <- K
  strata[["y"]] <- pheno
  strata[["Z"]] <- strat
  strata[["info"]] <- strat_cases
  strata[["ids"]] <- ids
  strata[["strat_miss"]] <- cases_nostrat
  strata[["sparse"]] <- sparse

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
  
  # Convert to percentiles (this becomes the "order" variable)
  percentiles <- ranks / length(ranks)
  
  # Add minimal jitter to break remaining ties for group assignment
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
  
  return(list(
    groups = as.integer(groups),
    order = percentiles  # Return the percentile ranks as the order
  ))
}

#' Create list of stratified phenotypes
#'
#' @param pheno Full phenotype data.frame
#' @param strat_cases Stratification info for cases
#' @param K Number of strata
#' @return List of phenotype data.frames, one per stratum
create_strata_list <- function(pheno, strat_cases, K, cov = NULL, ids) {

  strata <- vector("list", K)
  names(strata) <- paste0("group", 1:K)

  # Prepare covariate data if provided
  if (!is.null(cov)) {
    cov_df <- cov[match(ids, cov[, 1]), -(1:2), drop = FALSE]
    cov_df[] <- lapply(cov_df, function(x) {
      x <- as.numeric(x)
      x[is.na(x)] <- mean(x, na.rm = TRUE)
      x
    })
    rownames(cov_df) <- ids
  }

  for (k in 1:K) {
    # Select controls (pheno == 0) and cases from stratum k
    cases_in_stratum <- strat_cases[strat_cases$groups == k, 1]
    stratum_pheno <- pheno[pheno[, 3] == 0 | pheno[, 1] %in% cases_in_stratum, ]

    # Regress covariates from phenotype if provided
    if (!is.null(cov)) {
      y <- stratum_pheno[, 3]
      names(y) <- stratum_pheno[, 1]

      # Match covariate rows to stratum samples
      stratum_ids <- stratum_pheno[, 1]
      stratum_cov <- cov_df[match(stratum_ids, ids), , drop = FALSE]

      fit <- lm(y ~ ., data = stratum_cov)
      stratum_pheno[match(names(residuals(fit)), stratum_pheno[, 1]), 3] <- residuals(fit)
    }

    # Normalize phenotype (mean 0, sd 1)
    stratum_pheno[, 3] <- as.numeric(scale(stratum_pheno[, 3]))
    strata[[k]] <- stratum_pheno
  }

  return(strata)
}