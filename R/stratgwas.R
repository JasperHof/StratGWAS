#' Compute transformed phenotype using StratGWAS
#'
#' Compute genetic covariance of subgroups using Haseman-Elston regression
#'
#' @useDynLib StratGWAS, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' 
#' @param pheno Matrix containing binary outcome in PLINK format
#' @param strat_cont Matrix containing continuous stratification variables in PLINK format
#' @param strat_cat Matrix containing categorical stratification variables in PLINK format
#' @param cov Matrix containing covariates in PLINK format
#' @param filename Prefix of input .bed file
#' @param block_size Block size for HE regression
#' @param cor_g Pre-computed genetic covariance matrix (optional)
#' @param alpha Offset for stratification variables
#' @param max_criterion Maximum allowed inflation criterion
#' @param greedy_selection Use greedy variable selection (default TRUE)
#' @param min_improvement Minimum relative improvement for variable selection (default 0.05 = 5%)
#' @param verbose Print detailed progress information
#' @return Returns list with transformed phenotype and diagnostic information
#' @export
stratgwas <- function(pheno, filename, strat_cont = NULL, strat_cat = NULL, cov = NULL,
                      block_size = 500, cor_g = NULL, alpha = 0,
                      max_criterion = 0.01, greedy_selection = TRUE,
                      min_improvement = 0.05, verbose = TRUE) {
  # Validate input data
  validate_stratgwas_inputs(pheno, filename, strat_cont, strat_cat,
                            cov, block_size, cor_g, alpha)

  fam <- read.table(paste0(filename, ".fam"))

  # reduce phenotype
  pheno <- pheno[!is.na(pheno[, 3]), ]
  pheno <- pheno[pheno[, 1] %in% fam[, 1], ]
  ids <- pheno[, 1]

  # construct multivariate phenotype for HE-regression
  multi <- pheno[, 3]

  # Track variable names for reporting
  var_names <- c("Binary_Trait")

  # add continuous variables
  if(!is.null(strat_cont)){
    # reducing continuous stratification variable
    strat_cont <- strat_cont[strat_cont[, 1] %in% pheno[which(pheno[, 3] == 1), 1], ]

    # Get all continuous columns (excluding first 2 ID columns)
    cont_cols <- strat_cont[, -(1:2), drop = FALSE]

    # Initialize list to store spline basis for each continuous variable
    spline_list <- list()

    # Process each continuous column
    for(col_idx in 1:ncol(cont_cols)){
      col_data <- cont_cols[, col_idx]

      # Handle missing values - impute with mean
      col_data <- as.numeric(col_data)
      col_data[is.na(col_data)] <- mean(col_data, na.rm = TRUE)

      # Standardize
      col_data <- as.numeric(scale(col_data))

      # Apply spline transformation
      spline_mat <- spline_manual(col_data)
      spline_list[[col_idx]] <- spline_mat

      # Add variable names for each spline basis
      var_names <- c(var_names, paste0("Cont", col_idx, "_spline", 1:ncol(spline_mat)))
    }

    # Combine all spline matrices
    if(length(spline_list) > 0){
      all_splines <- do.call(cbind, spline_list)
      n_splines <- ncol(all_splines)
      
      # Add spline columns to multi
      multi <- cbind(multi, matrix(0, nrow = nrow(pheno), ncol = n_splines))
      multi[match(strat_cont[, 1], ids), (ncol(multi) - n_splines + 1):ncol(multi)] <- alpha + all_splines
    }
  }

  # add categorical variables
  if(!is.null(strat_cat)){
    # reducing continuous stratification variable
    strat_cat <- strat_cat[strat_cat[, 1] %in% pheno[which(pheno[, 3] == 1), 1], ]

    cat_cols <- strat_cat[, -(1:2), drop = FALSE]
    dummy_list <- list()

    # Process each categorical column
    for(col_idx in 1:ncol(cat_cols)){
      col_data <- cat_cols[, col_idx]

      # Handle missing values - replace with most common category (not ideal)
      if(any(is.na(col_data))){
        most_common <- names(sort(table(col_data), decreasing = TRUE))[1]
        col_data[is.na(col_data)] <- most_common
      }

      # Convert to factor if not already
      col_data <- as.factor(col_data)

      # Get unique levels (categories)
      levels_vec <- levels(col_data)
      n_levels <- length(levels_vec)

      # Create dummy variables
      if(n_levels > 1){
        dummy_mat <- matrix(0, nrow = length(col_data), ncol = n_levels - 1)
        for(i in 2:n_levels){
          dummy_mat[, i - 1] <- as.numeric(col_data == levels_vec[i])
        }
        dummy_list[[col_idx]] <- dummy_mat

        # Add variable names for each dummy variable
        var_names <- c(var_names, paste0("Cat", col_idx, "_", levels_vec[-1]))
      }
    }

    # Combine all dummy matrices
    if(length(dummy_list) > 0){
      all_dummies <- do.call(cbind, dummy_list)
      n_dummies <- ncol(all_dummies)
      
      # Add dummy columns to multi
      multi <- cbind(multi, matrix(0, nrow = nrow(pheno), ncol = n_dummies))
      multi[match(strat_cat[, 1], ids), (ncol(multi) - n_dummies + 1):ncol(multi)] <- alpha + all_dummies
    }
  }

  rownames(multi) <- ids
  colnames(multi) <- var_names

  # regress covariates from phenotype
  if(!is.null(cov)){
    cov_df <- cov[match(ids, cov[,1]), -(1:2), drop = F]
    cov_df[] <- lapply(cov_df, function(x) {
        x <- as.numeric(x)
        x[is.na(x)] <- mean(x, na.rm = TRUE)
        x
    })
    rownames(cov_df) <- ids

    # regress from all columns of multi
    for(i in 1:ncol(multi)){
      y <- multi[, i]
      names(y) <- ids

      fit <- lm(y ~ ., data = cov_df)
      multi[match(names(residuals(fit)), ids), i] <- residuals(fit)
    }
  }

  # normalize all columns
  for(i in 1:dim(multi)[2]) multi[, i] <- as.numeric(scale(multi[, i]))

  # compute genetic covariance using randomized HE regression
  if(is.null(cor_g)) cor_g <- he_multi_part(filename, multi, block_size)

  ############
  # continue computing correlation and transformation variable
  ############

  rg <- cor_g

  # check for negative heritability
  hers <- diag(cor_g)
  her_neg <- which(hers < 0)

  for(k in 1:dim(rg)[1]){
    if(hers[k] < 0) {
      cat(sprintf("SNP heritability of trait %d is negative, so will not compute genetic correlation \n", k))
      rg[k, ] <- NA
      rg[, k] <- NA
    } else {
      rg[k, ] <- rg[k, ] / sqrt(hers[k])
      rg[, k] <- rg[, k] / sqrt(hers[k])
    }
  }

  # store genetic correlation between binary trait and auxiliary variable
  if(ncol(multi) >= 2) {
    gamma <- rg[1, 2]
    h2_Z <- cor_g[2, 2]
  } else {
    gamma <- NA
    h2_Z <- NA
  }

  # compute environmental correlation
  cor_total = cor(multi, use = "complete.obs")
  cor_e = cor_total - cor_g

  if(length(her_neg) > 0){
    message("Heritability of term(s) ", paste(her_neg, collapse = ", "),
            " is negative and will be ignored.")

    if(length(her_neg) >= 3){
      warning("Many terms have negative heritability.")
    }

    # remove negative heritability terms from analysis
    cor_g_temp <- cor_g[-her_neg, -her_neg]
    cor_e_temp <- cor_e[-her_neg, -her_neg]
    multi_temp <- multi[, -her_neg]
    var_names_temp <- var_names[-her_neg]
  } else {
    # use full matrices if no negative heritabilities
    cor_g_temp <- cor_g
    cor_e_temp <- cor_e
    multi_temp <- multi
    var_names_temp <- var_names
  }

  # perform double decomposition
  #P <- eigen(cor_e_use)
  #S_inv <- P$vectors %*% diag(1/sqrt(P$values)) %*% t(P$vectors)
  #S <- P$vectors %*% diag(sqrt(P$values)) %*% t(P$vectors)
  #M <- S_inv %*% as.matrix(cor_g_use) %*% S_inv
  #V <- eigen(M)$vectors
  #U1 <- t(V) %*% S_inv

  #get weights
  #weights <- U1[1,]

  # compute estimated new heritability
  #D <- U1 %*% t(cor_g_use) %*% t(U1)
  #cat(sprintf("SNP heritability of transformed phenotype: %.4f\n\n", diag(D)[1]))

  # Apply greedy variable selection if requested and we have multiple variables
  if(greedy_selection && ncol(multi_temp) >= 2) {
    
    if(verbose) {
      cat("\n=== Greedy Variable Selection ===\n\n")
    }
    
    selection_result <- select_vars_greedy(
      cor_g = cor_g_temp,
      cor_e = cor_e_temp,
      min_improvement = min_improvement,
      verbose = verbose
    )
    
    # Use only selected variables
    selected_idx <- selection_result$selected_vars
    cor_g_use <- cor_g_temp[selected_idx, selected_idx, drop = FALSE]
    cor_e_use <- cor_e_temp[selected_idx, selected_idx, drop = FALSE]
    multi_use <- multi_temp[, selected_idx, drop = FALSE]
    var_names_use <- var_names_temp[selected_idx]
    
    # Get weights (for selected variables only)
    weights_selected <- selection_result$weights_selected
    
    # Create full weight vector (with zeros for unselected variables)
    weights_full <- numeric(ncol(multi_temp))
    weights_full[selected_idx] <- weights_selected
    
    h2_transformed <- selection_result$h2_final
    
    if(verbose) {
      cat(sprintf("\nSNP heritability of transformed phenotype: %.4f\n", h2_transformed))
    }
    
  } else {
    # No variable selection - use all variables
    cor_g_use <- cor_g_temp
    cor_e_use <- cor_e_temp
    multi_use <- multi_temp
    var_names_use <- var_names_temp
    selected_idx <- 1:ncol(multi_temp)
    
    if(ncol(multi_temp) == 1) {
      # Only binary trait
      weights_selected <- 1
      weights_full <- 1
      h2_transformed <- cor_g_use[1, 1]
      
    } else {
      # perform double decomposition
      P <- eigen(cor_e_use)
      S_inv <- P$vectors %*% diag(1/sqrt(P$values)) %*% t(P$vectors)
      S <- P$vectors %*% diag(sqrt(P$values)) %*% t(P$vectors)
      M <- S_inv %*% as.matrix(cor_g_use) %*% S_inv
      V <- eigen(M)$vectors
      U1 <- t(V) %*% S_inv

      # get weights
      weights_selected <- U1[1,]
      weights_full <- weights_selected

      # compute estimated new heritability
      D <- U1 %*% t(cor_g_use) %*% t(U1)
      h2_transformed <- as.numeric(diag(D)[1])
    }
    
    if(verbose) {
      cat(sprintf("SNP heritability of transformed phenotype: %.4f\n\n", h2_transformed))
    }
  }

# Compute genetic correlation between transformed and original trait
  if(ncol(multi_use) >= 2) {
    rg_result <- compute_rg_transformed(cor_g_use, weights_selected, binary_idx = 1)
    
    if(verbose) {
      cat(sprintf("Genetic correlation (transformed vs. original): %.4f\n", rg_result$rg))
      
      # Compute expected power improvement
      power_improvement <- compute_power_improvement(
        h2_original = rg_result$h2_binary,
        h2_transformed = rg_result$h2_transformed,
        rg = rg_result$rg
      )
      
      cat(sprintf("Expected GWAS power improvement: %.2fx (%.1f%% gain)\n\n",
                  power_improvement$improvement_factor,
                  power_improvement$expected_power_gain))
    }
    
    # Warning if genetic correlation is too low
    if(abs(rg_result$rg) < 0.7) {
      warning("Genetic correlation between transformed and original trait is low (", 
              round(rg_result$rg, 3), "). The transformed trait may not capture the ",
              "same genetic signal.")
    }
    
    rg_transformed <- rg_result$rg
    power_factor <- power_improvement$improvement_factor
    
  } else {
    rg_transformed <- 1.0  # Perfect correlation with itself
    power_factor <- 1.0
  }

  # compute transformed phenotype
  trans_pheno <- cbind(ids, ids, multi_use %*% weights_selected)

  # compute inflation criterion
  if(is.na(gamma) || ncol(multi) < 2){
    if(ncol(multi) < 2) {
      message("Cannot compute inflation criterion with only binary trait.")
    } else {
      message("Cannot compute inflation criterion due to negative heritability.")
    }
    criterion <- NA
  } else {
    # compute a2 - the correlation between stratification Y' and Z
    y_trans <- as.numeric(scale(trans_pheno[,3]))
    y <- as.numeric(scale(multi[,1]))
    z <- as.numeric(scale(multi[,2]))

    # fit regression
    fit <- lm(y_trans ~ y + z - 1)
    coefs <- fit$coefficients
    a2 <- coefs["z"]

    # Compute inflation criterion using ORIGINAL genetic correlation
    rg_original <- cor_g
    hers_original <- diag(cor_g)
    
    for(k in 1:nrow(rg_original)) {
      if(hers_original[k] > 0) {
        rg_original[k, ] <- rg_original[k, ] / sqrt(hers_original[k])
        rg_original[, k] <- rg_original[, k] / sqrt(hers_original[k])
      }
    }
    
    gamma_original <- rg_original[1, 2]
    h2_Z_original <- cor_g[2, 2]

    criterion = a2^2 * (1 - gamma_original^2) * h2_Z_original
    message(paste0("Expected inflation criterion is ", round(criterion, 4)))
    if(criterion > max_criterion) {
      warning(sprintf("Inflation criterion (%.4f) is greater than %.2f.", 
                     criterion, max_criterion))
    }
  }

  # Return list with information
  object = vector("list")
  object[["pheno"]] <- trans_pheno
  object[["multi"]] <- multi
  object[["cor_g"]] <- cor_g             # Original genetic correlation
  object[["cor_g_scaled"]] <- cor_g_use  # Used genetic correlation
  object[["cor_e"]] <- cor_e
  object[["rg"]] <- rg
  object[["criterion"]] <- criterion
  object[["weights"]] <- weights_full           # Weights for all variables (zeros for unselected)
  object[["weights_selected"]] <- weights_selected  # Weights for selected variables only

  return(object)
}