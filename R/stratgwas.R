#' Compute transformed phenotype using StratGWAS
#'
#' Compute genetic covariance of subgroups using Haseman-Elston regression
#'
#' @useDynLib StratGWAS, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' 
#' @param pheno Matrix containing binary outcome in PLINK format
#' @param strat Matrix containing stratification variable in PLINK format
#' @param filename Prefix of input .bed file
#' @return Returns covariance matrix of the strata
#' @export
stratgwas <- function(pheno, filename, strat_cont = NULL, strat_cat = NULL, cov = NULL,
                      block_size = 500, cor_g = NULL, alpha = 0,
                      constrain_inflation = TRUE, max_criterion = 0.01) {
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
  gamma <- rg[1, 2]
  h2_Z <- cor_g[2, 2]

  # compute environmental correlation
  cor_total = cor(multi, use = "complete.obs")
  cor_e = cor_total - cor_g

  if(length(her_neg) > 0){
    message("Heritability of term(s) ", paste(her_neg, collapse = ", "),
            " is negative and will be ignored.")

    if(length(her_neg) >= 3){
      warning("Too many terms have negative heritability.")
    }

    # remove negative heritability terms from analysis
    cor_g_use <- cor_g[-her_neg, -her_neg]
    cor_e_use <- cor_e[-her_neg, -her_neg]
    multi_use <- multi[, -her_neg]
  } else {
    # use full matrices if no negative heritabilities
    cor_g_use <- cor_g
    cor_e_use <- cor_e
    multi_use <- multi
  }

  # perform double decomposition
  P <- eigen(cor_e_use)
  S_inv <- P$vectors %*% diag(1/sqrt(P$values)) %*% t(P$vectors)
  S <- P$vectors %*% diag(sqrt(P$values)) %*% t(P$vectors)
  M <- S_inv %*% as.matrix(cor_g_use) %*% S_inv
  V <- eigen(M)$vectors
  U1 <- t(V) %*% S_inv

  # compute estimated new heritability
  D <- U1 %*% t(cor_g_use) %*% t(U1)
  cat(sprintf("SNP heritability of transformed phenotype: %.4f\n\n", diag(D)[1]))

  # compute transformed phenotype
  trans_pheno <- cbind(ids, ids, multi_use %*% weights)

  # compute inflation criterion
  if(is.na(gamma)){
    message("Can not compute inflation criterion due to negative heritability.")
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
  object[["cor_g"]] <- cor_g  # Original genetic correlation
  object[["cor_g_scaled"]] <- cor_g_use  # Scaled genetic correlation (if constrained)
  object[["cor_e"]] <- cor_e
  object[["rg"]] <- rg
  object[["criterion"]] <- criterion
  object[["weights"]] <- weights
  object[["constrained"]] <- constrain_inflation
  object[["scaling_applied"]] <- scaling_applied

  return(object)
}