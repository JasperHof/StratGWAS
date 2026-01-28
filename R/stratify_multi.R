#' Stratify cases with multiple variables
#'
#' Compute case subgroups by stratification based on multiple input variables
#'
#' @param pheno Binary input phenotype in PLINK format
#' @param strat_cont Continuous stratification variables in PLINK format
#' @param strat_bin Binary/categorical stratification variables in PLINK format
#' @param K Number of groups for continuous variables (default 5)
#' @param cov Covariates in PLINK format
#' @return Returns list containing subgroup phenotypes
#' @export
stratify_multi <- function(pheno, strat_cont = NULL, strat_bin = NULL, K = 5, cov = NULL) {

  # Store IDs
  ids <- pheno[, 1]
  control_ids <- pheno[which(pheno[, 3] == 0), 1]

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
  
  # Initialize for tracking all strata
  all_strat_info <- list()
  strata <- list()
  K_total <- 0
  all_cases_nostrat <- c()
  
  # Process continuous stratification variables
  if(!is.null(strat_cont)) {
    n_cont <- ncol(strat_cont) - 2
    
    for(i in 1:n_cont) {
      # Extract column
      strat_col <- strat_cont[, 2 + i]
      
      # Identify cases with missing stratification
      cases_nostrat <- strat_cont[which(strat_cont[, 1] %in% cases & is.na(strat_col)), 1]
      all_cases_nostrat <- c(all_cases_nostrat, cases_nostrat)
      
      # Extract stratification variable for cases only
      strat_cases <- strat_cont[which(pheno[, 3] == 1 & !(pheno[, 1] %in% cases_nostrat)), c(1, 2, 2 + i)]
      colnames(strat_cases) <- c("FID", "IID", "strat_val")
      strat_cases_vals <- strat_cases[, 3]
      
      # Determine stratification approach
      n_unique <- length(unique(strat_cases_vals[!is.na(strat_cases_vals)]))
      sparse <- n_unique < 10
      
      # Define groups
      if (sparse) {
        message(sprintf("Continuous variable %d: Only %d unique values detected. Using K=%d strata",
                        i, n_unique, n_unique))
        K_use <- n_unique
        strat_cases$groups <- match(strat_cases_vals, sort(unique(strat_cases_vals)))
        strat_cases$order <- NA
      } else {
        K_use <- K
        result <- assign_to_quantiles(strat_cases_vals, K)
        strat_cases$groups <- result$groups
        strat_cases$order <- result$order
      }
      
      # Create stratified phenotype lists for this variable
      for(k in 1:K_use) {
        K_total <- K_total + 1
        cases_in_stratum <- strat_cases[strat_cases$groups == k, 1]
        stratum_pheno <- pheno[pheno[, 3] == 0 | pheno[, 1] %in% cases_in_stratum, ]
        
        # Regress covariates if provided
        if(!is.null(cov)) {
          stratum_pheno <- apply_covariate_regression_single(stratum_pheno, cov, ids)
        }
        
        # Store in strata list
        strata[[paste0("group", K_total)]] <- stratum_pheno
      }
      
      # Store stratification info
      all_strat_info[[length(all_strat_info) + 1]] <- list(
        data = strat_cases,
        type = "continuous",
        K = K_use,
        cases_nostrat = cases_nostrat,
        var_index = i,
        group_start = K_total - K_use + 1,
        group_end = K_total
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
      all_cases_nostrat <- c(all_cases_nostrat, cases_nostrat)
      
      # Extract stratification variable for cases only
      strat_cases <- strat_bin[which(pheno[, 3] == 1 & !(pheno[, 1] %in% cases_nostrat)), c(1, 2, 2 + i)]
      colnames(strat_cases) <- c("FID", "IID", "strat_val")
      strat_cases_vals <- strat_cases[, 3]
      
      # Get unique categories
      categories <- sort(unique(strat_cases_vals[!is.na(strat_cases_vals)]))
      C <- length(categories)
      
      message(sprintf("Categorical variable %d: %d categories detected", i, C))
      
      # Assign groups based on categories
      strat_cases$groups <- match(strat_cases_vals, categories)
      strat_cases$order <- NA
      
      # Create stratified phenotype lists for this variable
      for(c in 1:C) {
        K_total <- K_total + 1
        cases_in_stratum <- strat_cases[strat_cases$groups == c, 1]
        stratum_pheno <- pheno[pheno[, 3] == 0 | pheno[, 1] %in% cases_in_stratum, ]
        
        # Regress covariates if provided
        if(!is.null(cov)) {
          stratum_pheno <- apply_covariate_regression_single(stratum_pheno, cov, ids)
        }
        
        # Store in strata list
        strata[[paste0("group", K_total)]] <- stratum_pheno
      }
      
      # Store stratification info
      all_strat_info[[length(all_strat_info) + 1]] <- list(
        data = strat_cases,
        type = "categorical",
        K = C,
        categories = categories,
        cases_nostrat = cases_nostrat,
        var_index = i,
        group_start = K_total - C + 1,
        group_end = K_total
      )
    }
  }
  
  # Check that we have at least one stratification variable
  if(length(all_strat_info) == 0) {
    stop("No stratification variables provided. Must provide strat_cont or strat_bin.")
  }
  
  # Create a combined info dataframe for backward compatibility
  # This combines all stratification info, but individuals can appear multiple times
  all_cases_combined <- do.call(rbind, lapply(all_strat_info, function(x) {
    df <- x$data
    # Adjust group numbers to be global
    df$groups <- df$groups + x$group_start - 1
    df
  }))
  
  # Create a multivariate dataframe for linear regression + SumHer
  multi <- cbind(ids, ids)

  for(k in 1:length(all_strat_info)){
    var_add  <- NULL

    for(j in all_strat_info[[k]]$group_start:all_strat_info[[k]]$group_end){
      add <- rep(NA, length(ids))
      add[which(ids %in% control_ids)] <- 0
      add[which(ids %in% all_cases_combined[all_cases_combined$groups == j, 1])] <- 1
      var_add <- cbind(var_add, add)
    }

    colnames(var_add) <- paste0(all_strat_info[[k]]$type, "_", all_strat_info[[k]]$var_index, "_", 1:all_strat_info[[k]]$K)

    # In the case of a continuous variable, add this
    if(all_strat_info[[k]]$type == "continuous"){
      add_cont <- rep(NA, length(ids))
      add_cont[match(all_strat_info[[k]]$data[, 1], ids)] <- all_strat_info[[k]]$data[, 3]
      var_add <- cbind(var_add, add_cont)

      colnames(var_add)[ncol(var_add)] <- paste0(all_strat_info[[k]]$type, "_", all_strat_info[[k]]$var_index)
    }

    multi <- cbind(multi, var_add)
  }

  # Determine if covariates were used
  cov_used <- !is.null(cov)

  # Return list with information (matching original structure)
  strata[["K"]] <- K_total
  strata[["y"]] <- pheno
  strata[["multi"]] <- multi
  strata[["info"]] <- all_cases_combined
  strata[["ids"]] <- ids
  strata[["strat_miss"]] <- all_cases_nostrat
  strata[["cov_used"]] <- cov_used
  strata[["strat_details"]] <- all_strat_info  # Detailed info about each variable
  strata[["n_cont"]] <- if(!is.null(strat_cont)) ncol(strat_cont) - 2 else 0
  strata[["n_bin"]] <- if(!is.null(strat_bin)) ncol(strat_bin) - 2 else 0
  
  return(strata)
}

#' Apply covariate regression to a single phenotype
#'
#' @param stratum_pheno Phenotype dataframe for one stratum
#' @param cov Covariate matrix
#' @param ids All individual IDs
#' @return Phenotype dataframe with residualized values
apply_covariate_regression_single <- function(stratum_pheno, cov, ids) {
  cov_df <- cov[match(ids, cov[, 1]), -(1:2), drop = FALSE]
  cov_df[] <- lapply(cov_df, function(x) {
    x <- as.numeric(x)
    x[is.na(x)] <- mean(x, na.rm = TRUE)
    x
  })
  rownames(cov_df) <- ids
  
  y <- stratum_pheno[, 3]
  names(y) <- stratum_pheno[, 1]
  
  # Match covariate rows to stratum samples
  stratum_ids <- stratum_pheno[, 1]
  stratum_cov <- cov_df[match(stratum_ids, ids), , drop = FALSE]
  
  fit <- lm(y ~ ., data = stratum_cov)
  stratum_pheno[match(names(residuals(fit)), stratum_pheno[, 1]), 3] <- residuals(fit)
  
  return(stratum_pheno)
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