#' Transform phenotype based on genetic covariance
#'
#' Create a transformed phenotype that maximizes power to detect genetic
#' associations by leveraging the genetic covariance structure between strata.
#' This is the third step in the StratGWAS workflow.
#'
#' @useDynLib StratGWAS, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' 
#' @param strata Object returned from \code{\link{stratify}}
#' @param gencov Genetic covariance object returned by \code{\link{compute_gencov}}
#' @param outfile Name/path prefix for output files
#'   
#' @return Returns a list containing:
#'   \item{transformed_pheno}{Data frame with FID, IID, and transformed phenotype values}
#'   \item{weights}{Named vector of transformation weights for each stratum}
#'   \item{inflation_criteria}{Data frame with expected inflation statistics for 
#'     continuous variables, including a2 (regression coefficient), genetic correlation 
#'     with binary phenotype, heritability of stratification variable, and expected 
#'     inflation factor}
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' data(pheno)
#' data(strat_cont)
#' filename_bed <- system.file("extdata", "data.bed", package = "StratGWAS")
#' filename <- gsub(".bed", "", filename_bed)
#' outfile <- tempfile("transform")
#'
#' strata <- stratify(pheno, strat_cont = strat_cont, K = 5)
#' gencov <- compute_gencov(strata, filename, nr_blocks = 1000, outfile)
#' trans <- transform(strata, gencov)
#'
#' # Examine transformation results
#' head(trans$transformed_pheno)      # Transformed phenotype values
#' print(trans$weights)               # Weights for each stratum
#' print(trans$inflation_criteria)    # Expected inflation (continuous vars only)
#'
#' # Check distribution of transformed phenotype
#' hist(trans$transformed_pheno[, 3], main = "Transformed Phenotype Distribution")
#'
#' # With categorical variables
#' data(strat_cat)
#' strata_cat <- stratify(pheno, strat_cat = strat_cat)
#' gencov_cat <- compute_gencov(strata_cat, filename, nr_blocks = 1000, outfile)
#' trans_cat <- transform(strata_cat, gencov_cat)
#' }
#' 
#' @export
transform <- function(strata, gencov, outfile, spar = 0.8, smooth = TRUE) {
  
  # Get ids
  ids <- strata$y[, 2]
  control_ids <- strata$y[which(strata$y[, 3] == 0), 2]

  # Check smoothing parameter
  if (strata$K < 4) {
    message("Need K > 3 to enable smoothing - will continue without smoothing.")
    smooth <- FALSE
  }

  # Information about prevalence
  total_prev <- mean(as.numeric(strata$y[, 3]), na.rm = T)
  cont_prev <- (total_prev / strata$K) / (1 - total_prev + total_prev / strata$K)

  # Genetic covariance
  gencov_all <- gencov$gencov

  # Extract strata and compute eigendecomposition
  names <- c()
  for(k in 1:length(strata$strat_details)) names <- c(names, paste0(strata$strat_details[[k]]$type, "_", strata$strat_details[[k]]$var_index,"_", 1:strata$strat_details[[k]]$K))

  idx <- which(colnames(gencov_all) %in% names)
  gencov_use <- gencov_all[idx, idx]
  multi_use <- strata$multi[, colnames(strata$multi) %in% names]

  # Update observed scale heritabilities of categorical groups
  for(k in 1:length(strata$strat_details)){
    if(strata$strat_details[[k]]$type == "categorical"){
      vars <- paste0(strata$strat_details[[k]]$type,"_",strata$strat_details[[k]]$var_index,"_",1:strata$strat_details[[k]]$K)
      h2_obs <- diag(gencov_use[vars, vars])

      # Ensure we have a numeric matrix
      multi_subset <- as.matrix(strata$multi[, vars, drop = FALSE])
      multi_subset <- apply(multi_subset, 2, as.numeric)
      prevs <- colMeans(multi_subset, na.rm = TRUE)

      ### scale heritability estimates
      z_from <- dnorm(qnorm(1 - prevs))
      z_to <- dnorm(qnorm(1 - cont_prev))

      scale_from <- (prevs * (1 - prevs)) / z_from^2
      scale_to <- (cont_prev * (1 - cont_prev)) / z_to^2

      fac <- scale_from
      #fac <- scale_from / scale_to
      ###

      # compute h2 relative to prevalence of strata for continuous variables
      #fac <- (cont_prev * (1 - cont_prev)) / (prevs * (1 - prevs))
      h2_adj <- h2_obs * fac

      # update genetic covariance
      for(i in 1:length(vars)){
        gencov_use[vars[i], ] <- gencov_use[vars[i], ] * sqrt(fac[i])
        gencov_use[, vars[i]] <- gencov_use[, vars[i]] * sqrt(fac[i])
      }
    } else {
      vars <- paste0(strata$strat_details[[k]]$type, "_", strata$strat_details[[k]]$var_index,"_", 1:strata$strat_details[[k]]$K)
      h2_obs <- diag(gencov_use[vars, vars])

      # Ensure we have a numeric matrix
      multi_subset <- as.matrix(strata$multi[, vars, drop = FALSE])
      multi_subset <- apply(multi_subset, 2, as.numeric)
      prevs <- colMeans(multi_subset, na.rm = TRUE)

      # Scale heritability estimates
      z_from <- dnorm(qnorm(1 - prevs))
      scale <- (prevs * (1 - prevs)) / z_from^2
      h2_adj <- h2_obs * scale

      # update genetic covariance
      for(i in 1:length(vars)){
        gencov_use[vars[i], ] <- gencov_use[vars[i], ] * sqrt(scale[i])
        gencov_use[, vars[i]] <- gencov_use[, vars[i]] * sqrt(scale[i])
      }
    }
  }

  # Smooth matrix if needed? Turn off for now
  #if(min(eigen(gencov_use)$values) <= 0) {
  #  gencov_use <- as.matrix(Matrix::nearPD(gencov_use, corr = FALSE)$mat)
  #}

  # Compute eigenvector transformation
  trans <- eigen(gencov_use)$vectors[, 1]
  trans <- trans / sqrt(sum(trans^2))  # normalize
  if (mean(trans) < 0) trans <- -trans
  names(trans) <- names

  # get jack-knife SE estimates of weights
  weights_all <- NULL
  weights_means <- NULL

  if (gencov$jack) {
    weights_se <- weights_se_jack(strata, gencov, trans, gencov_all, names)
    weights_uni <- weights_univariate(strata, gencov, trans, gencov_all, names)
  } else {
    weights_se <- NULL
    weights_uni <- NULL
  }

  # Initialize transformed phenotype with controls
  trans_pheno <- data.frame(FID = ids, IID = ids, Pheno = 0)

  # Compute transformed phenotype separately for continuous and categorical
  for(k in 1:length(strata$strat_details)){

    # Get weights for this stratification variable
    weights <- trans[names(trans) %in% paste0(strata$strat_details[[k]]$type,"_",strata$strat_details[[k]]$var_index,"_",1:strata$strat_details[[k]]$K)]
    weights_all <- c(weights_all, weights)

    if (strata$strat_details[[k]]$type == "continuous" && smooth) {
      
      # Compute medians for (non-)scaled stratification variable
      strat_scale <- as.numeric(scale(strata$strat_details[[k]]$data[, 3]))
      strat_noscale <- as.numeric(strata$strat_details[[k]]$data[, 3])
      
      medians_obs <- unlist(lapply(1:strata$strat_details[[k]]$K, function(x) mean(strat_scale[which(strata$strat_details[[k]]$data$groups == x)], na.rm = TRUE)))
      medians_obs_noscale <- unlist(lapply(1:strata$strat_details[[k]]$K, function(x) mean(strat_noscale[which(strata$strat_details[[k]]$data$groups == x)], na.rm = TRUE)))
      weights_means <- c(weights_means, medians_obs_noscale)

      # Smooth through values
      fit <- smooth.spline(medians_obs, as.numeric(weights), spar = spar)
      trans_pred <- predict(fit, strat_scale)$y

      # Create a variable to add to transformation variable
      trans_add <- rep(NA, nrow(trans_pheno))
      trans_add[ids %in% control_ids] <- 0
      trans_add[match(strata$strat_details[[k]]$data[, 2], ids)] <- trans_pred
      trans_add[is.na(trans_add)] <- mean(trans_pred, na.rm = T) # individuals with missing value for the stratificatin variable

    } else {

      # Create subset of multi matrix for this variable
      multi_use_sub <- multi_use[, colnames(multi_use) %in% paste0(strata$strat_details[[k]]$type,"_",strata$strat_details[[k]]$var_index,"_",1:strata$strat_details[[k]]$K), drop = FALSE]
      multi_use_sub <- apply(multi_use_sub, 2, as.numeric)
      miss <- which(rowSums(is.na(multi_use_sub)) == ncol(multi_use_sub))

      # NA values for non-continuous stratification variables
      weights_means <- c(weights_means, rep(NA, ncol(multi_use_sub)))

      multi_use_sub[is.na(multi_use_sub)] <- 0
      trans_add <- multi_use_sub %*% weights
      
      # Impute missing phenotypes by mean
      if (length(miss) > 0) {
        trans_add[miss] <- mean(trans_add[rowSums(multi_use_sub == 0) != ncol(multi_use_sub)])
      }
    }

    trans_pheno[, 3] <- trans_pheno[, 3] + trans_add

  }

  # Write trans_pheno and weights to output file
  if (gencov$jack == TRUE) {
    weights_df <- data.frame("weights" = weights_all, "weights_SE" = weights_se, 
                             "weights_uni" = weights_uni$weights_uni, "weights_uni_SE" = weights_uni$weights_uni_se)
  } else {
    weights_df <- data.frame("weights" = weights_all)
  }
  weights_df$cont_means <- weights_means

  message(paste0("Writing transformed phenotype to ",
                 paste0(outfile, ".transformed")))
  write.table(trans_pheno, paste0(outfile, ".transformed"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  message(paste0("Writing StratGWAS weights to ",
                 paste0(outfile, ".weights")))
  write.table(weights_df, paste0(outfile, ".weights"),
              quote = FALSE, row.names = TRUE, col.names = TRUE)

  # Compute the inflation factors - need a2, gencor, and h2_Z
  h2_y <- gencov_all[1, 1]

  inflation_results <- data.frame(
    variable = character(),
    a2 = numeric(),
    rg = numeric(),
    h2_Z = numeric(),
    expected_inflation = numeric(),
    stringsAsFactors = FALSE
  )

  for(k in 1:length(strata$strat_details)){
    if(strata$strat_details[[k]]$type == "continuous"){

      var_name <- paste0(strata$strat_details[[k]]$type, "_", strata$strat_details[[k]]$var_index)

      idx <- which(colnames(gencov_all) == var_name)
      h2_Z <- gencov_all[idx, idx]

      if(h2_y * h2_Z < 0){
        message(paste0("Can not compute expected inflation criterion of ", var_name,
        " due to negative h2 estimate of binary trait or stratification variable"))
      }

      # Prepare variables
      Z <- strata$strat_details[[k]]$original[, 3]
      Z[is.na(Z)] <- mean(Z, na.rm = TRUE)
      
      y <- as.numeric(strata$y[, 3])
      y[is.na(y)] <- mean(y, na.rm = TRUE)

      # Scale variables
      Z <- as.numeric(scale(Z))
      y <- as.numeric(scale(y))
      Y_trans <- as.numeric(scale(trans_pheno[, 3]))

      # Fit regression to compute a2
      fit <- lm(Y_trans ~ y + Z - 1)
      coefs <- fit$coefficients
      a2 <- coefs["Z"]

      # Compute inflation criterion using genetic correlation
      rg <- gencov_all[1, idx] / sqrt(h2_y * h2_Z)
      exp_inflation <- a2^2 * (1 - rg^2) * h2_Z

      message(sprintf("Inflation criterion of %s is %.4f",
                      var_name, exp_inflation))

      inflation_results <- rbind(
        inflation_results,
        data.frame(
          variable = var_name,
          a2 = a2,
          rg = rg,
          h2_Z = h2_Z,
          expected_inflation = exp_inflation,
          stringsAsFactors = FALSE
        )
      )
    }
  }

  cat("\n")

  return(list(
    transformed_pheno = trans_pheno,
    weights = trans,
    weights_se = weights_se,
    weights_uni = weights_uni$weights_uni,
    weights_uni_se = weights_uni$weights_uni_se,
    inflation_criteria = inflation_results
  ))
}
