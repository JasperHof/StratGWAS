#' Stratify cases
#'
#' Compute case subgroups by stratification based on input variable
#'
#' @param phenofile Binary input phenotype in PLINK format
#' @param stratfile Stratification variable in PLINK format
#' @param filename Prefix of input .bed file
#' @return Returns list containing subgroup phenotypes
#' @export
stratify_multi <- function(pheno, strat_cont = NULL, strat_bin = NULL, K = 5, cov = NULL) {

  # Check input data
  #stratify_checks_multi(pheno, strat, K)

  # Store IDs
  ids <- pheno[, 1]

  # Convert to data.frame
  pheno <- as.data.frame(pheno)

  if(!is.null(strat_cont)){
    strat_cont <- as.data.frame(strat_cont)
    strat_cont <- strat_cont[match(ids, strat_cont[, 1]), ]
    strat_cont[, 1] <- strat_cont[, 2] <- ids
  }
  if(!is.null(strat_bin)){
    strat_bin <- as.data.frame(strat_bin)
    strat_bin <- strat_bin[match(ids, strat_bin[, 1]), ]
    strat_bin[, 1] <- strat_bin[, 2] <- ids
  }

  # Identify cases
  cases <- pheno[which(pheno[, 3] == 1), 1]
  
  # Initialize list to store all stratified phenotypes
  all_strata <- list()
  strat_info_list <- list()
  K_list <- c()
  var_names <- c()
  
  # Process continuous stratification variables
  if(!is.null(strat_cont)) {
    n_cont <- ncol(strat_cont) - 2
    
    for(i in 1:n_cont) {
      # Extract column
      strat_col <- strat_cont[, 2 + i]
      
      # Identify cases with missing stratification
      cases_nostrat <- strat_cont[which(strat_cont[, 1] %in% cases & is.na(strat_col)), 1]
      
      # Extract stratification variable for cases only
      strat_cases <- strat_cont[which(pheno[, 3] == 1 & !(pheno[, 1] %in% cases_nostrat)), ]
      strat_cases_vals <- strat_cases[, 2 + i]
      
      # Determine stratification approach
      n_unique <- length(unique(strat_cases_vals[!is.na(strat_cases_vals)]))
      sparse <- n_unique < 10
      
      # Define groups
      if (sparse) {
        # Groups defined by unique values
        message(sprintf("Continuous variable %d: Only %d unique values detected. Using K=%d strata",
                        i, n_unique, n_unique))
        K_use <- n_unique
        strat_cases$groups <- match(strat_cases_vals, sort(unique(strat_cases_vals)))
      } else {
        # Quantile-based stratification
        K_use <- K
        result <- assign_to_quantiles(strat_cases_vals, K)
        strat_cases$groups <- result$groups
        strat_cases$order <- result$order
      }
      
      # Create stratified phenotype lists for this variable
      strata_i <- create_strata_list(pheno, strat_cases, K_use, cov = cov, ids = ids)
      
      # Store with informative names
      for(k in 1:K_use) {
        var_name <- paste0("cont", i, "_group", k)
        all_strata[[var_name]] <- strata_i[[k]]
        var_names <- c(var_names, var_name)
      }
      
      # Store metadata
      K_list <- c(K_list, K_use)
      strat_info_list[[paste0("cont", i)]] <- list(
        type = "continuous",
        strat_cases = strat_cases,
        K = K_use,
        sparse = sparse,
        cases_nostrat = cases_nostrat
      )
    }
  }
  
  # Process binary/categorical stratification variables
  if(!is.null(strat_bin)) {
    n_bin <- ncol(strat_bin) - 2
    
    for(i in 1:n_bin) {
      # Extract column
      strat_col <- strat_bin[, 2 + i]
      
      # Identify cases with missing stratification
      cases_nostrat <- strat_bin[which(strat_bin[, 1] %in% cases & is.na(strat_col)), 1]
      
      # Extract stratification variable for cases only
      strat_cases <- strat_bin[which(pheno[, 3] == 1 & !(pheno[, 1] %in% cases_nostrat)), ]
      strat_cases_vals <- strat_cases[, 2 + i]
      
      # Get unique categories
      categories <- sort(unique(strat_cases_vals[!is.na(strat_cases_vals)]))
      C <- length(categories)
      
      message(sprintf("Categorical variable %d: %d categories detected", i, C))
      
      # Assign groups based on categories
      strat_cases$groups <- match(strat_cases_vals, categories)
      
      # Create stratified phenotype lists for this variable
      strata_i <- create_strata_list(pheno, strat_cases, C, cov = cov, ids = ids)
      
      # Store with informative names
      for(c in 1:C) {
        var_name <- paste0("bin", i, "_cat", c)
        all_strata[[var_name]] <- strata_i[[c]]
        var_names <- c(var_names, var_name)
      }
      
      # Store metadata
      K_list <- c(K_list, C)
      strat_info_list[[paste0("bin", i)]] <- list(
        type = "categorical",
        strat_cases = strat_cases,
        K = C,
        categories = categories,
        cases_nostrat = cases_nostrat
      )
    }
  }
  
  # Check that we have at least one stratification variable
  if(length(all_strata) == 0) {
    stop("No stratification variables provided. Must provide strat_cont or strat_bin.")
  }
  
  # Determine if covariates were used
  cov_used <- !is.null(cov)
  
  # Return list with information
  result <- list()
  result[["strata"]] <- all_strata
  result[["var_names"]] <- var_names
  result[["K_list"]] <- K_list
  result[["strat_info"]] <- strat_info_list
  result[["y"]] <- pheno
  result[["strat_cont"]] <- strat_cont
  result[["strat_bin"]] <- strat_bin
  result[["ids"]] <- ids
  result[["cov_used"]] <- cov_used
  result[["n_cont"]] <- if(!is.null(strat_cont)) ncol(strat_cont) - 2 else 0
  result[["n_bin"]] <- if(!is.null(strat_bin)) ncol(strat_bin) - 2 else 0
  result[["total_strata"]] <- length(all_strata)
  
  return(result)
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