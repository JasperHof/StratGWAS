#' @keywords internal  
validate_stratgwas_inputs <- function(pheno, filename, strat_cont = NULL, 
                                      strat_cat = NULL, cov = NULL,
                                      block_size = 500, cor_g = NULL, 
                                      alpha = 0) {
  
  # Check pheno
  if(is.null(pheno)) {
    stop("pheno cannot be NULL")
  }
  
  if(!is.matrix(pheno) && !is.data.frame(pheno)) {
    stop("pheno must be a matrix or data.frame")
  }
  
  if(ncol(pheno) < 3) {
    stop("pheno must have at least 3 columns: FID, IID, and phenotype")
  }
  
  # Check phenotype is binary (0/1 or 1/2)
  pheno_vals <- unique(pheno[, 3])
  pheno_vals <- pheno_vals[!is.na(pheno_vals)]
  
  if(!all(pheno_vals %in% c(0, 1))) {
    stop("Phenotype column (column 3) must contain binary values (0/1)")
  }
  
  if(sum(pheno[, 3] == 1, na.rm = TRUE) < 10) {
    stop("Too few cases in phenotype. Need at least 10 cases")
  }
  
  # Check for duplicate IDs
  if(any(duplicated(pheno[, 1]))) {
    stop("Duplicate IDs found in pheno")
  }
  
  # Check filename and associated files
  if(is.null(filename) || !is.character(filename)) {
    stop("filename must be a character string")
  }
  
  required_files <- c(paste0(filename, ".bed"),
                     paste0(filename, ".bim"),
                     paste0(filename, ".fam"))
  
  missing_files <- required_files[!file.exists(required_files)]
  if(length(missing_files) > 0) {
    stop("Missing required files: ", paste(missing_files, collapse = ", "))
  }
  
  # Check strat_cont
  if(!is.null(strat_cont)) {
    if(!is.matrix(strat_cont) && !is.data.frame(strat_cont)) {
      stop("strat_cont must be a matrix or data.frame")
    }
    
    if(ncol(strat_cont) < 3) {
      stop("strat_cont must have at least 3 columns: FID, IID, and at least one continuous variable")
    }
    
    # Check that continuous columns are numeric
    cont_cols <- strat_cont[, -(1:2), drop = FALSE]
    for(i in 1:ncol(cont_cols)) {
      if(!is.numeric(as.numeric(cont_cols[, i]))) {
        stop("Column ", i + 2, " in strat_cont cannot be coerced to numeric")
      }
    }
    
    # Check for overlap with pheno IDs
    overlap <- sum(strat_cont[, 1] %in% pheno[, 1])
    if(overlap == 0) {
      stop("No overlapping IDs between strat_cont and pheno")
    }
    if(overlap < nrow(strat_cont) * 0.5) {
      warning("Less than 50% of strat_cont IDs found in pheno")
    }
  }
  
  # Check strat_cat
  if(!is.null(strat_cat)) {
    if(!is.matrix(strat_cat) && !is.data.frame(strat_cat)) {
      stop("strat_cat must be a matrix or data.frame")
    }
    
    if(ncol(strat_cat) < 3) {
      stop("strat_cat must have at least 3 columns: FID, IID, and at least one categorical variable")
    }
    
    # Check that categorical columns have reasonable number of levels
    cat_cols <- strat_cat[, -(1:2), drop = FALSE]
    for(i in 1:ncol(cat_cols)) {
      n_unique <- length(unique(cat_cols[!is.na(cat_cols[, i]), i]))
      
      if(n_unique == 0) {
        stop("Column ", i + 2, " in strat_cat has no non-missing values")
      }
      
      if(n_unique == 1) {
        warning("Column ", i + 2, " in strat_cat has only one unique value")
      }
      
      if(n_unique > 10) {
        warning("Column ", i + 2, " in strat_cat has ", n_unique,
                " levels. Are you sure this is categorical?")
      }
    }
    
    # Check for overlap with pheno IDs
    overlap <- sum(strat_cat[, 1] %in% pheno[, 1])
    if(overlap == 0) {
      stop("No overlapping IDs between strat_cat and pheno")
    }
    if(overlap < nrow(strat_cat) * 0.5) {
      warning("Less than 50% of strat_cat IDs found in pheno")
    }
  }
  
  # Check that at least one stratification variable is provided
  if(is.null(strat_cont) && is.null(strat_cat)) {
    stop("Must provide at least one of strat_cont or strat_cat")
  }
  
  # Check cov
  if(!is.null(cov)) {
    if(!is.matrix(cov) && !is.data.frame(cov)) {
      stop("cov must be a matrix or data.frame")
    }
    
    if(ncol(cov) < 3) {
      stop("cov must have at least 3 columns: FID, IID, and at least one covariate")
    }
    
    # Check for overlap with pheno IDs
    overlap <- sum(cov[, 1] %in% pheno[, 1])
    if(overlap == 0) {
      stop("No overlapping IDs between cov and pheno")
    }
    if(overlap < nrow(cov) * 0.5) {
      warning("Less than 50% of cov IDs found in pheno")
    }
  }
  
  # Check block_size
  if(!is.numeric(block_size) || block_size <= 50) {
    stop("block_size must be a positive number > 50")
  }
  
  if(block_size < 100) {
    warning("block_size is very small (< 100). This may be inefficient.")
  }
  
  # Check cor_g
  if(!is.null(cor_g)) {
    if(!is.matrix(cor_g)) {
      stop("cor_g must be a matrix")
    }
    
    if(nrow(cor_g) != ncol(cor_g)) {
      stop("cor_g must be a square matrix")
    }
    
    if(!isSymmetric(cor_g)) {
      warning("cor_g is not symmetric")
    }
  }
  
  # Check alpha
  if(!is.numeric(alpha) || length(alpha) != 1) {
    stop("alpha must be a single numeric value")
  }
  
  if(abs(alpha) > 10) {
    warning("alpha value is large (|alpha| > 10).")
  }
  
  # All checks passed
  invisible(NULL)
}

#' @keywords internal  
stratify_checks <- function(pheno, strat, K) {
  
  # Check 1: pheno validation
  if (missing(pheno) || is.null(pheno)) {
    stop("'pheno' argument is required and cannot be NULL")
  }
  
  if (!is.data.frame(pheno) && !is.matrix(pheno)) {
    stop("'pheno' must be a data.frame or matrix in PLINK format (FID, IID, phenotype)")
  }
  
  if (ncol(pheno) < 3) {
    stop("'pheno' must have at least 3 columns (FID, IID, phenotype)")
  }
  
  if (nrow(pheno) == 0) {
    stop("'pheno' has no rows")
  }
  
  # Check phenotype values (should be binary: 0 = control, 1 = case, optionally -9 = missing)
  pheno_vals <- unique(pheno[, 3])
  pheno_vals <- pheno_vals[!is.na(pheno_vals)]
  
  if (!all(pheno_vals %in% c(0, 1))) {
    warning("'pheno' column 3 contains values other than 0/1. ")
  }
  
  # Check for at least some cases and controls
  n_cases <- sum(pheno[, 3] == 1, na.rm = TRUE)
  n_controls <- sum(pheno[, 3] == 0, na.rm = TRUE)
  n_missing <- sum(is.na(pheno[, 3]) | pheno[, 3] == -9 | pheno[, 3] == 2)
  
  if (n_cases == 0) {
    stop("'pheno' contains no cases (phenotype = 1)")
  }
  
  if (n_controls == 0) {
    stop("'pheno' contains no controls (phenotype = 0)")
  }
  
  if (n_missing > 0) {
    message("Note: 'pheno' contains ", n_missing, " individuals with missing phenotype")
  }
  
  # Check 2: strat validation
  if (missing(strat) || is.null(strat)) {
    stop("'strat' argument is required and cannot be NULL")
  }
  
  if (!is.data.frame(strat) && !is.matrix(strat)) {
    stop("'strat' must be a data.frame or matrix in PLINK format (FID, IID, stratification_value)")
  }
  
  if (ncol(strat) < 3) {
    stop("'strat' must have at least 3 columns (FID, IID, stratification_value)")
  }
  
  if (nrow(strat) == 0) {
    stop("'strat' has no rows")
  }
  
  # Check for duplicate IDs in pheno
  if (any(duplicated(pheno[, 1]))) {
    stop("'pheno' contains duplicate IDs in column 1 (FID)")
  }
  
  # Check for duplicate IDs in strat
  if (any(duplicated(strat[, 1]))) {
    stop("'strat' contains duplicate IDs in column 1 (FID)")
  }

  # Check for number of unique values in stratification variable
  cases <- pheno[which(pheno[,3] == 1),1]
  strat_vals <- table(strat[which(strat[,1] %in% cases),3])
  
  if(max(strat_vals)/sum(strat_vals) > 0.10){
    warning("Some stratification variables are identical for >10% of cases")
  }
  if(length(strat_vals)/sum(strat_vals) < 0.10){
    warning("Low number of identical values for stratification variable")
  }

  if (any(duplicated(strat[, 1]))) {
    stop("'strat' contains duplicate IDs in column 1 (FID)")
  }

  # Check 3: K validation
  if (missing(K)) {
    stop("'K' argument is required")
  }
  
  if (!is.numeric(K) || length(K) != 1) {
    stop("'K' must be a single numeric value")
  }
  
  if (K < 2 || K != round(K)) {
    stop("'K' must be an integer >= 2 (number of strata)")
  }
  
  if (K > 20) {
    warning("'K' = ", K, " is quite large. Are you sure you want this many strata?")
  }
  
  # Check 4: Overlap between pheno and strat
  pheno_ids_cases <- pheno[which(pheno[,3] == 1), 1]
  pheno_ids <- pheno[, 1]
  strat_ids <- strat[, 1]
  
  ids_in_both <- intersect(pheno_ids, strat_ids)
  ids_only_pheno <- setdiff(pheno_ids_cases, strat_ids)
  ids_only_strat <- setdiff(strat_ids, pheno_ids)
  
  if (length(ids_in_both) == 0) {
    stop("No overlapping IDs between 'pheno' and 'strat'. Check that ID columns match.")
  }
  
  if (length(ids_only_strat) > 0) {
    message("Note: ", length(ids_only_strat), " IDs in 'strat' are not in 'pheno' (will be ignored)")
  }
  
  # Check 5: Stratification variable among cases
  # Match strat to pheno to check cases
  strat_matched <- strat[match(pheno[, 1], strat[, 1]), ]
  case_mask <- which(pheno[, 3] == 1)
  strat_cases <- strat_matched[case_mask, 3]
  
  n_cases_with_strat <- sum(!is.na(strat_cases))
  n_cases_missing_strat <- sum(is.na(strat_cases))
  
  if (n_cases_with_strat == 0) {
    stop("No cases have non-missing stratification values. Cannot stratify.")
  }
  
  if (n_cases_missing_strat > 0) {
    warning(n_cases_missing_strat, " cases (", 
            round(100 * n_cases_missing_strat / n_cases, 1), 
            "%) have missing stratification values. ",
            "These will be imputed with the median.")
  }
  
  # Check that there are enough cases per stratum
  if (n_cases_with_strat < K) {
    stop("Not enough cases with stratification values (", n_cases_with_strat, 
         ") to create K = ", K, " strata")
  }
  
  min_cases_per_stratum <- n_cases_with_strat / K
  if (min_cases_per_stratum < 10) {
    warning("With ", n_cases_with_strat, " cases divided into ", K, " strata, ",
            "you'll have ~", round(min_cases_per_stratum, 1), " cases per stratum on average. ",
            "Consider using fewer strata (smaller K) for more stable estimates.")
  }
  
  # Check 6: Stratification variable properties
  strat_vals <- strat_cases[!is.na(strat_cases)]
  
  if (!is.numeric(strat_vals)) {
    stop("Stratification variable (column 3 of 'strat') must be numeric")
  }
  
  # Check for constant stratification variable
  if (length(unique(strat_vals)) == 1) {
    stop("Stratification variable has no variation (all values are ", strat_vals[1], ")")
  }
  
  # Check for very limited variation
  if (length(unique(strat_vals)) < K) {
    stop("Stratification variable has only ", length(unique(strat_vals)), 
         " unique values among cases, but K = ", K, " strata requested. ",
         "Reduce K or use a stratification variable with more variation.")
  }
  
  # Warn if distribution looks problematic
  if (length(unique(strat_vals)) < 2 * K) {
    warning("Stratification variable has only ", length(unique(strat_vals)), 
            " unique values among cases for K = ", K, " strata. ",
            "This may result in uneven stratum sizes.")
  }
  
  return(TRUE)
}

#' @keywords internal
compute_gencov_checks <- function(strata, filename, nr_blocks, outfile, 
                                   SumHer, ss_list, lds) {
  
  # Check 1: strata object validation
  if (missing(strata) || is.null(strata)) {
    stop("'strata' argument is required and cannot be NULL")
  }
  
  if (!is.list(strata)) {
    stop("'strata' must be a list object returned from stratify()")
  }
  
  required_strata_fields <- c("K", "y", "info")
  missing_fields <- setdiff(required_strata_fields, names(strata))
  if (length(missing_fields) > 0) {
    stop("'strata' object is missing required fields: ", 
         paste(missing_fields, collapse = ", "))
  }
  
  # Check K is valid
  if (!is.numeric(strata$K) || strata$K < 3 || strata$K != round(strata$K)) {
    stop("'strata$K' must be an integer >= 3")
  }
  
  K <- strata$K
  
  # Check that group1, group2, ..., groupK exist
  for (k in 1:K) {
    group_name <- paste0("group", k)
    if (!group_name %in% names(strata)) {
      stop("'strata' is missing '", group_name, "' (expected K=", K, " groups)")
    }
    
    # Check group structure
    if (!is.data.frame(strata[[group_name]]) && !is.matrix(strata[[group_name]])) {
      stop("'strata$", group_name, "' must be a data.frame or matrix")
    }
    
    if (ncol(strata[[group_name]]) < 3) {
      stop("'strata$", group_name, "' must have at least 3 columns (FID, IID, phenotype)")
    }
  }

  # Check 2: Genotype files exist
  if (missing(filename) || is.null(filename) || !is.character(filename)) {
    stop("'filename' must be a character string specifying the .bed file prefix")
  }

  required_files <- paste0(filename, c(".bed", ".bim", ".fam"))
  missing_files <- required_files[!file.exists(required_files)]
  if (length(missing_files) > 0) {
    stop("Required genotype files not found: ", paste(missing_files, collapse = ", "))
  }

  # Check 3: nr_blocks validation
  if (!is.numeric(nr_blocks) || nr_blocks < 1 || nr_blocks != round(nr_blocks)) {
    stop("'nr_blocks' must be a positive integer")
  }

  # Check 4: outfile validation
  if (missing(outfile) || is.null(outfile) || !is.character(outfile)) {
    stop("'outfile' must be a character string")
  }

  # Check if output directory exists
  outdir <- dirname(outfile)
  if (outdir != "." && !dir.exists(outdir)) {
    stop("Output directory does not exist: ", outdir)
  }

  # Check 5: SumHer parameter
  if (!is.logical(SumHer) || length(SumHer) != 1) {
    stop("'SumHer' must be a single logical value (TRUE or FALSE)")
  }

  # Check 6: ss_list validation (if provided)
  if (!is.null(ss_list)) {
    if (!is.list(ss_list)) {
      stop("'ss_list' must be a list")
    }
    
    K_tot <- K + 2
    if (length(ss_list) != K_tot) {
      stop("'ss_list' must have length K+2 = ", K_tot, ", but has length ", 
           length(ss_list))
    }
    
    # Check each summary statistic data frame
    required_ss_cols <- c("Predictor", "N", "Chisq", "Beta", "SE")
    for (i in 1:K_tot) {
      if (!is.data.frame(ss_list[[i]])) {
        stop("'ss_list[[", i, "]]' must be a data.frame")
      }
      
      missing_cols <- setdiff(required_ss_cols, names(ss_list[[i]]))
      if (length(missing_cols) > 0) {
        stop("'ss_list[[", i, "]]' is missing columns: ", 
             paste(missing_cols, collapse = ", "))
      }
      
      # Check for NAs in critical columns
      if (any(is.na(ss_list[[i]]$N))) {
        stop("'ss_list[[", i, "]]$N' contains NA values")
      }
      
      # Check all ss have same number of SNPs
      if (i > 1 && nrow(ss_list[[i]]) != nrow(ss_list[[1]])) {
        stop("All summary statistics must have the same number of SNPs. ",
             "ss_list[[1]] has ", nrow(ss_list[[1]]), " rows, but ",
             "ss_list[[", i, "]] has ", nrow(ss_list[[i]]), " rows")
      }
    }
  }
  
  # Check 7: lds validation (if provided)
  if (!is.null(lds)) {
    if (!is.data.frame(lds)) {
      stop("'lds' must be a data.frame")
    }
    
    required_ld_cols <- c("Predictor", "Tagging")
    missing_cols <- setdiff(required_ld_cols, names(lds))
    if (length(missing_cols) > 0) {
      stop("'lds' is missing required columns: ", 
           paste(missing_cols, collapse = ", "))
    }
    
    if (any(is.na(lds$Tagging))) {
      stop("'lds$Tagging' contains NA values")
    }
    
    if (any(lds$Tagging < 0)) {
      stop("'lds$Tagging' contains negative values")
    }
  }
  
  return(TRUE)
}

#' @keywords internal
compute_gencov_checks <- function(strata, filename, nr_blocks, outfile, 
                                   SumHer, ss_list, lds) {
  
  # Check 1: strata object validation
  if (missing(strata) || is.null(strata)) {
    stop("'strata' argument is required and cannot be NULL")
  }
  
  if (!is.list(strata)) {
    stop("'strata' must be a list object returned from stratify()")
  }
  
  required_strata_fields <- c("K", "y", "info")
  missing_fields <- setdiff(required_strata_fields, names(strata))
  if (length(missing_fields) > 0) {
    stop("'strata' object is missing required fields: ", 
         paste(missing_fields, collapse = ", "))
  }
  
  # Check K is valid
  if (!is.numeric(strata$K) || strata$K < 3 || strata$K != round(strata$K)) {
    stop("'strata$K' must be an integer >= 3")
  }
  
  K <- strata$K
  
  # Check that group1, group2, ..., groupK exist
  for (k in 1:K) {
    group_name <- paste0("group", k)
    if (!group_name %in% names(strata)) {
      stop("'strata' is missing '", group_name, "' (expected K=", K, " groups)")
    }
    
    # Check group structure
    if (!is.data.frame(strata[[group_name]]) && !is.matrix(strata[[group_name]])) {
      stop("'strata$", group_name, "' must be a data.frame or matrix")
    }
    
    if (ncol(strata[[group_name]]) < 3) {
      stop("'strata$", group_name, "' must have at least 3 columns (FID, IID, phenotype)")
    }
  }

  # Check 2: Genotype files exist
  if (missing(filename) || is.null(filename) || !is.character(filename)) {
    stop("'filename' must be a character string specifying the .bed file prefix")
  }

  required_files <- paste0(filename, c(".bed", ".bim", ".fam"))
  missing_files <- required_files[!file.exists(required_files)]
  if (length(missing_files) > 0) {
    stop("Required genotype files not found: ", paste(missing_files, collapse = ", "))
  }

  # Check 3: nr_blocks validation
  if (!is.numeric(nr_blocks) || nr_blocks < 1 || nr_blocks != round(nr_blocks)) {
    stop("'nr_blocks' must be a positive integer")
  }

  # Check 4: outfile validation
  if (missing(outfile) || is.null(outfile) || !is.character(outfile)) {
    stop("'outfile' must be a character string")
  }

  # Check if output directory exists
  outdir <- dirname(outfile)
  if (outdir != "." && !dir.exists(outdir)) {
    stop("Output directory does not exist: ", outdir)
  }

  # Check 5: SumHer parameter
  if (!is.logical(SumHer) || length(SumHer) != 1) {
    stop("'SumHer' must be a single logical value (TRUE or FALSE)")
  }

  # Check 6: ss_list validation (if provided)
  if (!is.null(ss_list)) {
    if (!is.list(ss_list)) {
      stop("'ss_list' must be a list")
    }
    
    K_tot <- K + 2
    if (length(ss_list) != K_tot) {
      stop("'ss_list' must have length K+2 = ", K_tot, ", but has length ", 
           length(ss_list))
    }
    
    # Check each summary statistic data frame
    required_ss_cols <- c("Predictor", "N", "Chisq", "Beta", "SE")
    for (i in 1:K_tot) {
      if (!is.data.frame(ss_list[[i]])) {
        stop("'ss_list[[", i, "]]' must be a data.frame")
      }
      
      missing_cols <- setdiff(required_ss_cols, names(ss_list[[i]]))
      if (length(missing_cols) > 0) {
        stop("'ss_list[[", i, "]]' is missing columns: ", 
             paste(missing_cols, collapse = ", "))
      }
      
      # Check for NAs in critical columns
      if (any(is.na(ss_list[[i]]$N))) {
        stop("'ss_list[[", i, "]]$N' contains NA values")
      }
      
      # Check all ss have same number of SNPs
      if (i > 1 && nrow(ss_list[[i]]) != nrow(ss_list[[1]])) {
        stop("All summary statistics must have the same number of SNPs. ",
             "ss_list[[1]] has ", nrow(ss_list[[1]]), " rows, but ",
             "ss_list[[", i, "]] has ", nrow(ss_list[[i]]), " rows")
      }
    }
  }
  
  # Check 7: lds validation (if provided)
  if (!is.null(lds)) {
    if (!is.data.frame(lds)) {
      stop("'lds' must be a data.frame")
    }
    
    required_ld_cols <- c("Predictor", "Tagging")
    missing_cols <- setdiff(required_ld_cols, names(lds))
    if (length(missing_cols) > 0) {
      stop("'lds' is missing required columns: ", 
           paste(missing_cols, collapse = ", "))
    }
    
    if (any(is.na(lds$Tagging))) {
      stop("'lds$Tagging' contains NA values")
    }
    
    if (any(lds$Tagging < 0)) {
      stop("'lds$Tagging' contains negative values")
    }
  }
  
  return(TRUE)
}

#' @keywords internal
spline_manual <- function(x){
  x <- as.numeric(x)

  # quantiles as interior knots
  k <- quantile(x, probs = c(1/3, 2/3))

  tp <- function(z) pmax(z, 0)^3

  B <- cbind(
    x,
    x^2,
    x^3,
    tp(x - k[1]),
    tp(x - k[2])
  )

  return(B)
}