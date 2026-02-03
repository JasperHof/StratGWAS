#' Stratify cases based on input variables
#'
#' Compute case subgroups by stratification based on multiple input variables.
#' This function creates strata of cases based on continuous and/or categorical
#' variables, which can then be used for genetic covariance estimation and
#' phenotype transformation in stratified GWAS analysis.
#'
#' @param pheno Binary input phenotype in PLINK format (data frame or matrix with
#'   columns: FID, IID, phenotype). Phenotype should be coded as 0 (controls) and
#'   1 (cases). Missing values coded as NA.
#' @param strat_cont Optional: continuous stratification variables in PLINK format
#'   (data frame or matrix with columns: FID, IID, variable1, variable2, ...). 
#'   Each variable will be divided into K quantile-based groups.
#' @param strat_cat Optional: binary/categorical stratification variables in PLINK 
#'   format (data frame or matrix with columns: FID, IID, variable1, variable2, ...).
#'   Categories are automatically detected from unique values.
#' @param K Number of quantile groups for continuous variables (default: 5).
#'   Categorical variables use their natural number of categories.
#' 
#' @return Returns a list containing:
#'   \item{K}{Number of groups used for continuous variables}
#'   \item{y}{Input phenotype data frame}
#'   \item{multi}{Multi-phenotype matrix for genetic covariance estimation, with 
#'     binary indicators for each stratum and continuous variables}
#'   \item{info}{Combined information about all strata with group assignments}
#'   \item{ids}{Vector of all individual IDs}
#'   \item{strat_miss}{Vector of case IDs with missing stratification values}
#'   \item{strat_details}{Detailed list for each stratification variable containing:
#'     \itemize{
#'       \item \code{data}: Stratification data for cases
#'       \item \code{original}: Original stratification variable data (for continuous)
#'       \item \code{type}: "continuous" or "categorical"
#'       \item \code{K}: Number of groups/categories
#'       \item \code{categories}: Category values (for categorical)
#'       \item \code{cases_nostrat}: Cases with missing values
#'       \item \code{var_index}: Variable index
#'       \item \code{group_start}, \code{group_end}: Global group indices
#'     }}
#'   \item{n_cont}{Number of continuous stratification variables}
#'   \item{n_bin}{Number of categorical stratification variables}
#' 
#' @examples
#' \dontrun{
#' # Load example data
#' data(pheno)
#' data(strat_cont)
#' data(strat_cat)
#' 
#' # Example 1: Stratify by single continuous variable (e.g., age at diagnosis)
#' strata_cont <- stratify(pheno, strat_cont = strat_cont, K = 5)
#' 
#' # View stratification summary
#' table(strata_cont$info$groups)  # Cases per stratum
#' length(strata_cont$strat_miss)  # Cases with missing values
#' 
#' # Example 2: Stratify by categorical variable (e.g., disease subtype)
#' strata_cat <- stratify(pheno, strat_cat = strat_cat)
#' 
#' # Check detected categories
#' strata_cat$strat_details[[1]]$categories
#' strata_cat$strat_details[[1]]$K  # Number of categories
#' 
#' # Example 3: Stratify by both continuous and categorical variables
#' strata_both <- stratify(pheno, 
#'                         strat_cont = strat_cont, 
#'                         strat_cat = strat_cat, 
#'                         K = 5)
#' 
#' # Check structure
#' strata_both$n_cont  # Number of continuous variables
#' strata_both$n_bin   # Number of categorical variables
#' 
#' # View multi-phenotype matrix (first few rows and columns)
#' head(strata_both$multi[, 1:6])
#' 
#' # Example 4: Use different number of quantiles
#' strata_10 <- stratify(pheno, strat_cont = strat_cont, K = 10)
#' 
#' # Example 5: Complete workflow
#' filename_bed <- system.file("extdata", "data.bed", package = "StratGWAS")
#' filename <- gsub(".bed", "", filename_bed)
#' outfile <- tempfile("gwas_output")
#' 
#' # Stratify
#' strata <- stratify(pheno, strat_cont = strat_cont, K = 5)
#' 
#' # Compute genetic covariance
#' gencov <- compute_gencov(strata, filename, nr_blocks = 1000, outfile)
#' 
#' # Transform phenotype
#' trans <- transform(strata, gencov)
#' 
#' # Run GWAS
#' results <- linear(trans, filename, outfile, nr_blocks = 1000)
#' head(results$results)
#' }
#' 
#' @export
stratify <- function(pheno, strat_cont = NULL, strat_cat = NULL, K = 5) {

  # Check input data
  check_stratify_inputs(pheno, strat_cont, strat_cat, K)

  # Store IDs
  ids <- pheno[, 2]
  control_ids <- pheno[which(pheno[, 3] == 0), 2]

  # Convert to data.frame
  pheno <- as.data.frame(pheno)

  if(!is.null(strat_cont)){
    strat_cont <- as.data.frame(strat_cont)
    strat_cont <- strat_cont[match(ids, strat_cont[, 2]), ]
    strat_cont[, 1] <- strat_cont[, 2] <- ids
  }
  if(!is.null(strat_cat)){
    strat_cat <- as.data.frame(strat_cat)
    strat_cat <- strat_cat[match(ids, strat_cat[, 2]), ]
    strat_cat[, 1] <- strat_cat[, 2] <- ids
  }

  # Identify cases
  cases <- pheno[which(pheno[, 3] == 1), 2]
  
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
      cases_nostrat <- strat_cont[which(strat_cont[, 2] %in% cases & is.na(strat_col)), 2, drop = TRUE]
      all_cases_nostrat <- c(all_cases_nostrat, cases_nostrat)
      
      # Extract stratification variable for cases only
      strat_cases <- strat_cont[which(pheno[, 3] == 1 & !(pheno[, 2] %in% cases_nostrat)), c(1, 2, 2 + i)]
      colnames(strat_cases) <- c("FID", "IID", "strat_val")
      strat_cases_vals <- strat_cases[, 3]
      
      # Define groups
      K_use <- K
      K_total <- K_total + K_use
      result <- assign_to_quantiles(strat_cases_vals, K)
      strat_cases$groups <- result$groups
      strat_cases$order <- result$order
      
      # Store stratification info
      all_strat_info[[length(all_strat_info) + 1]] <- list(
        data = strat_cases,
        original = strat_cont[, c(1, 2, 2 + i)],
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
  if(!is.null(strat_cat)) {
    n_bin <- ncol(strat_cat) - 2
    
    for(i in 1:n_bin) {
      # Extract column
      strat_col <- strat_cat[, 2 + i]
      
      # Identify cases with missing stratification
      cases_nostrat <- strat_cat[which(strat_cat[, 2] %in% cases & is.na(strat_col)), 2]
      all_cases_nostrat <- c(all_cases_nostrat, cases_nostrat)
      
      # Extract stratification variable for cases only
      strat_cases <- strat_cat[which(pheno[, 3] == 1 & !(pheno[, 2] %in% cases_nostrat)), c(1, 2, 2 + i)]
      colnames(strat_cases) <- c("FID", "IID", "strat_val")
      strat_cases_vals <- strat_cases[, 3]
      
      # Get unique categories
      categories <- sort(unique(strat_cases_vals[!is.na(strat_cases_vals)]))
      C <- length(categories)
      
      message(sprintf("Categorical variable %d: %d categories detected", i, C))
      
      # Assign groups based on categories
      K_total <- K_total + C
      strat_cases$groups <- match(strat_cases_vals, categories)
      strat_cases$order <- NA
      
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
    stop("No stratification variables provided. Must provide strat_cont or strat_cat")
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
  multi <- cbind(ids, ids, pheno[, 3])
  colnames(multi) <- c("FID", "IID", "pheno")

  for(k in 1:length(all_strat_info)){
    var_add  <- NULL

    for(j in all_strat_info[[k]]$group_start:all_strat_info[[k]]$group_end){
      add <- rep(NA, length(ids))
      add[which(ids %in% control_ids)] <- 0
      add[which(ids %in% all_cases_combined[all_cases_combined$groups == j, 2])] <- 1
      var_add <- cbind(var_add, add)
    }

    colnames(var_add) <- paste0(all_strat_info[[k]]$type, "_", all_strat_info[[k]]$var_index, "_", 1:all_strat_info[[k]]$K)

    # In the case of a continuous variable, add this (INCLUDING ALL CASES!)
    if(all_strat_info[[k]]$type == "continuous"){
      add_cont <- rep(NA, length(ids))
      add_cont <- all_strat_info[[k]]$original[match(ids, all_strat_info[[k]]$original[, 2]), 3]
      var_add <- cbind(var_add, add_cont)

      colnames(var_add)[ncol(var_add)] <- paste0(all_strat_info[[k]]$type, "_", all_strat_info[[k]]$var_index)
    }

    multi <- cbind(multi, var_add)
  }

  cat("\n")

  # Return list with information (matching original structure)
  strata[["K"]] <- K
  strata[["y"]] <- pheno
  strata[["multi"]] <- multi
  strata[["info"]] <- all_cases_combined
  strata[["ids"]] <- ids
  strata[["strat_miss"]] <- all_cases_nostrat
  strata[["strat_details"]] <- all_strat_info  # Detailed info about each variable
  strata[["n_cont"]] <- if(!is.null(strat_cont)) ncol(strat_cont) - 2 else 0
  strata[["n_bin"]] <- if(!is.null(strat_cat)) ncol(strat_cat) - 2 else 0
  
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