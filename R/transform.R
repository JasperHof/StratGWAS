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
transform <- function(strata, gencov, outfile) {
  
  # <performs some checks here> #

  # Get ids
  ids <- strata$y[, 2]
  control_ids <- strata$y[which(strata$y[, 3] == 0), 2]

  # Genetic covariance
  gencov <- gencov$gencov

  # Extract strata and compute eigendecomposition
  names <- c()
  for(k in 1:length(strata$strat_details)) names <- c(names, paste0(strata$strat_details[[k]]$type,"_",strata$strat_details[[k]]$var_index,"_",1:strata$strat_details[[k]]$K))

  idx <- which(colnames(gencov) %in% names)
  gencov_use <- gencov[idx, idx]
  multi_use <- strata$multi[, colnames(strata$multi) %in% names]

  # Update observed scale heritabilities of categorical groups
  for(k in 1:length(strata$strat_details)){
    if(strata$strat_details[[k]]$type == "categorical"){
      vars <- paste0(strata$strat_details[[k]]$type,"_",strata$strat_details[[k]]$var_index,"_",1:strata$strat_details[[k]]$K)

      h2_obs <- diag(gencov_use[vars, vars])
      
      # Ensure we have a numeric matrix
      multi_subset <- as.matrix(multi_use[, vars, drop = FALSE])
      multi_subset <- apply(multi_subset, 2, as.numeric)
      prevs <- colMeans(multi_subset, na.rm = TRUE)

      # compute relative to prevalence 1 / K (used for cont. variables)
      fac <- (1/strata$K * (1 - 1/strata$K)) / (prevs * (1 - prevs))
      h2_adj <- h2_obs * fac

      # update genetic covariance
      for(i in 1:length(vars)){
        gencov_use[vars[i], ] <- gencov_use[vars[i], ] * sqrt(fac[i])
        gencov_use[, vars[i]] <- gencov_use[, vars[i]] * sqrt(fac[i])
      }
    } 
  } 

  # Compute eigenvector transformation
  trans <- eigen(gencov_use)$vectors[, 1]
  names(trans) <- names

  # Initialize transformed phenotype with controls
  trans_pheno <- data.frame(FID = ids, IID = ids, Pheno = 0)

  # Compute transformed phenotype separately for continuous and categorical
  for(k in 1:length(strata$strat_details)){

    # Get weights for this stratification variable
    weights <- trans[names(trans) %in% paste0(strata$strat_details[[k]]$type,"_",strata$strat_details[[k]]$var_index,"_",1:strata$strat_details[[k]]$K)]

    if(strata$strat_details[[k]]$type == "continuous"){
      
      # Compute medians for scaled stratification variable
      strat_scale <- as.numeric(scale(strata$strat_details[[k]]$data[, 3]))
      medians_obs <- unlist(lapply(1:strata$strat_details[[k]]$K, function(x) mean(strat_scale[which(strata$strat_details[[k]]$data$groups == x)], na.rm = TRUE)))

      # Smooth through values
      fit <- smooth.spline(medians_obs, as.numeric(weights), spar = 0.2)
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

      multi_use_sub[is.na(multi_use_sub)] <- 0
      trans_add <- multi_use_sub %*% weights
      
      # Impute missing phenotypes by mean
      if(length(miss) > 0){
        trans_add[miss] <- mean(trans_add[rowSums(multi_use_sub == 0) != ncol(multi_use_sub)])
      }
    }

    trans_pheno[, 3] <- trans_pheno[, 3] + trans_add

  }

  # Write trans_pheno to output file
  message(paste0("Writing transformed phenotype to ",
                 paste0(outfile, ".transformed")))
  write.table(trans_pheno, paste0(outfile, ".transformed"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)

  # Compute the inflation factors - need a2, gencor, and h2_Z
  h2_y <- gencov[1, 1]

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

      idx <- which(colnames(gencov) == var_name)
      h2_Z <- gencov[idx, idx]

      if(h2_y * h2_Z < 0){
        message(paste0("Can not compute expected inflation criterion of ", var_name,
        " due to negative h2 estimate of binary trait or stratification variable"))
      }

      # Prepare variables
      Z <- strata$strat_details[[k]]$original[, 3]
      Z[is.na(Z)] <- mean(Z, na.rm = TRUE)
      
      y <- strata$y[, 3]
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
      rg <- gencov[1, idx] / sqrt(h2_y * h2_Z)
      exp_inflation <- a2^2 * (1 - rg^2) * h2_Z

      message(sprintf("Expected inflation criterion of %s is %.4f",
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
    inflation_criteria = inflation_results
  ))
}
