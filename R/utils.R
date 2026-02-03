#' Check inputs for stratify function
#'
#' Validates all input parameters for the stratify function
#'
#' @param pheno Binary input phenotype in PLINK format
#' @param strat_cont Continuous stratification variables in PLINK format
#' @param strat_cat Categorical stratification variables in PLINK format
#' @param K Number of groups for continuous variables
#' @param cov Covariates in PLINK format
#' @return NULL if all checks pass, stops with error message otherwise
#' @keywords internal
check_stratify_inputs <- function(pheno, strat_cont, strat_cat, K) {
  
  # Check pheno
  if (is.null(pheno)) {
    stop("pheno cannot be NULL")
  }
  
  if (!is.matrix(pheno) && !is.data.frame(pheno)) {
    stop("pheno must be a matrix or data.frame")
  }
  
  if (ncol(pheno) < 3) {
    stop("pheno must have at least 3 columns (FID, IID, phenotype)")
  }
  
  if (nrow(pheno) == 0) {
    stop("pheno must have at least one row")
  }
  
  # Check phenotype column is binary
  pheno_vals <- unique(pheno[!is.na(pheno[, 3]), 3])
  if (!all(pheno_vals %in% c(0, 1))) {
    stop("Phenotype column (column 3) must contain only 0 (control) and 1 (case) values")
  }
  
  # Check that we have at least one case and one control
  n_cases <- sum(pheno[, 3] == 1, na.rm = TRUE)
  n_controls <- sum(pheno[, 3] == 0, na.rm = TRUE)
  
  if (n_cases == 0) {
    stop("pheno must contain at least one case (phenotype = 1)")
  }
  
  if (n_controls == 0) {
    stop("pheno must contain at least one control (phenotype = 0)")
  }
  
  # Check for duplicate IDs in pheno
  if (any(duplicated(pheno[, 2]))) {
    stop("Duplicate IDs found in pheno (column 2)")
  }
  
  # Check K parameter
  if (!is.numeric(K) || length(K) != 1) {
    stop("K must be a single numeric value")
  }
  
  if (K < 2) {
    stop("K must be at least 2")
  }
  
  if (K != floor(K)) {
    stop("K must be an integer")
  }
  
  if (K > n_cases) {
    stop(sprintf("K (%d) cannot be greater than the number of cases (%d)", K, n_cases))
  }
  
  # Check that at least one stratification variable is provided
  if (is.null(strat_cont) && is.null(strat_cat)) {
    stop("At least one stratification variable must be provided (strat_cont or strat_cat)")
  }
  
  # Check strat_cont
  if (!is.null(strat_cont)) {
    if (!is.matrix(strat_cont) && !is.data.frame(strat_cont)) {
      stop("strat_cont must be a matrix or data.frame")
    }
    
    if (ncol(strat_cont) < 3) {
      stop("strat_cont must have at least 3 columns (FID, IID, and at least one stratification variable)")
    }
    
    if (nrow(strat_cont) == 0) {
      stop("strat_cont must have at least one row")
    }
    
    # Check for duplicate IDs in strat_cont
    if (any(duplicated(strat_cont[, 2]))) {
      stop("Duplicate IDs found in strat_cont (column 2)")
    }
    
    # Check that stratification columns are numeric
    for (i in 3:ncol(strat_cont)) {
      col_vals <- strat_cont[!is.na(strat_cont[, i]), i]
      if (length(col_vals) > 0 && !is.numeric(col_vals)) {
        tryCatch({
          test_numeric <- as.numeric(col_vals)
          if (any(is.na(test_numeric) & !is.na(col_vals))) {
            stop(sprintf("Column %d in strat_cont contains non-numeric values", i))
          }
        }, warning = function(w) {
          stop(sprintf("Column %d in strat_cont contains non-numeric values", i))
        })
      }
    }
    
    # Check that we have enough cases with non-missing stratification values
    for (i in 3:ncol(strat_cont)) {
      # Get case IDs from pheno
      case_ids <- pheno[pheno[, 3] == 1, 2]
      
      # Match cases to strat_cont and check for non-missing values
      strat_cont_matched <- strat_cont[match(case_ids, strat_cont[, 2]), ]
      cases_with_strat <- sum(!is.na(strat_cont_matched[, i]))
      
      if (cases_with_strat < K) {
        stop(sprintf("Stratification variable %d in strat_cont has only %d cases with non-missing values, but K = %d", 
                     i - 2, cases_with_strat, K))
      }
    }
  }
  
  # Check strat_cat
  if (!is.null(strat_cat)) {
    if (!is.matrix(strat_cat) && !is.data.frame(strat_cat)) {
      stop("strat_cat must be a matrix or data.frame")
    }
    
    if (ncol(strat_cat) < 3) {
      stop("strat_cat must have at least 3 columns (FID, IID, and at least one stratification variable)")
    }
    
    if (nrow(strat_cat) == 0) {
      stop("strat_cat must have at least one row")
    }
    
    # Check for duplicate IDs in strat_cat
    if (any(duplicated(strat_cat[, 2]))) {
      stop("Duplicate IDs found in strat_cat (column 2)")
    }
    
    # Check that each categorical variable has at least 2 categories among cases
    for (i in 3:ncol(strat_cat)) {
      case_ids <- pheno[pheno[, 3] == 1, 2]
      
      # Match cases to strat_cat and get non-missing values
      strat_cat_matched <- strat_cat[match(case_ids, strat_cat[, 2]), ]
      case_strat_vals <- strat_cat_matched[!is.na(strat_cat_matched[, i]), i]
      
      n_categories <- length(unique(case_strat_vals))
      
      if (n_categories < 2) {
        stop(sprintf("Stratification variable %d in strat_cat has only %d unique category/categories among cases (need at least 2)", 
                     i - 2, n_categories))
      }
    }
  }
  
  # All checks passed
  return(invisible(NULL))
}