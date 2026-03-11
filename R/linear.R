#' Linear regression GWAS on transformed phenotype
#'
#' Perform genome-wide association analysis using the transformed phenotype.
#' Covariates can be included and will be regressed out before analysis.
#' This is the fourth and final step in the StratGWAS workflow.
#'
#' @useDynLib StratGWAS, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#'
#' @param trans Object returned from \code{\link{transform}}, containing the
#'   transformed phenotype and associated weights
#' @param filename Prefix of genotype .bed file (without the .bed extension)
#' @param outfile Name/path prefix for output files
#' @param nr_blocks Block size for reading in genotype data (default: 1000)
#' @param cov Optional: data frame with covariates in PLINK format (FID, IID,
#'   and covariate columns). Missing covariate values will be imputed to their
#'   column means. Covariates will be regressed out of the transformed phenotype
#'   before GWAS analysis.
#'
#' @return GWAS results are written to <outfile>.trans.pheno1
#'
#' @examples
#' \dontrun{
#' # Basic usage with covariates
#' data(pheno)
#' data(strat_cont)
#' data(cov)
#' filename_bed <- system.file("extdata", "data.bed", package = "StratGWAS")
#' filename <- gsub(".bed", "", filename_bed)
#' outfile <- tempfile("gwas")
#'
#' strata <- stratify(pheno, strat_cont = strat_cont, K = 5)
#' gencov <- compute_gencov(strata, filename, nr_blocks = 1000, outfile)
#' trans <- transform(strata, gencov)
#' linear(trans, filename, outfile, nr_blocks = 1000, cov = cov)
#'
#' gwas_results <- read.table(paste0(outfile,".trans.pheno1"))
#'
#' # View top hits
#' top_hits <- gwas_results$gwas_results[order(gwas_results$gwas_results$Pvalue), ]
#' head(top_hits, 20)
#' }
#'
#' @export
linear <- function(trans, filename, outfile, nr_blocks = 1000, cov = NULL) {

  # Extract transformed phenotype
  trans_pheno <- trans$transformed_pheno

  # Read genotype IDs
  fam_ids <- as.character(read.table(paste0(filename, ".fam"))[, 2])

  # Match transformed phenotype to genotype file order
  pheno_matched <- rep(NA_real_, length(fam_ids))
  names(pheno_matched) <- fam_ids

  matched_idx <- match(trans_pheno[, 2], fam_ids)
  valid <- !is.na(matched_idx)
  pheno_matched[matched_idx[valid]] <- trans_pheno[valid, 3]

  # Regress out covariates if provided
  if (!is.null(cov)) {
    # Match covariates to genotype file order
    cov_ids <- as.character(cov[, 2])
    cov_data <- as.matrix(cov[, -(1:2), drop = FALSE])  # ensure matrix, not data.frame

    match_idx <- match(fam_ids, cov_ids)  # for each fam_id, find its row in cov

    cov_matched <- cov_data[match_idx, , drop = FALSE]   # reorder in one step
    cov_matched <- apply(cov_matched, 2, as.numeric)
    rownames(cov_matched) <- fam_ids

    # Impute missing covariate values with their means
    for (k in 1:ncol(cov_matched)) {
      missing_idx <- is.na(cov_matched[, k])
      if (any(missing_idx)) {
        mean_val <- mean(cov_matched[, k], na.rm = TRUE)
        cov_matched[missing_idx, k] <- mean_val
        cat(sprintf("Imputed %d missing values for covariate %s with mean (%.4f)\n", 
                    sum(missing_idx), colnames(cov_matched)[k], mean_val))
      }
    }

    # Regress out covariates from phenotype (only for non-missing phenotypes)
    non_missing <- !is.na(pheno_matched)

    if (sum(non_missing) > ncol(cov_matched)) {
      fit <- lm(pheno_matched[non_missing] ~ cov_matched[non_missing, , drop = FALSE])
      pheno_matched[non_missing] <- residuals(fit)
      
      cat(sprintf("Regressed out %d covariate(s) from transformed phenotype\n", 
                  ncol(cov_matched)))
    } else {
      warning("Not enough observations to regress out covariates")
    }
  }
  
  # Scale residualized phenotype
  non_missing <- !is.na(pheno_matched)
  if (sum(non_missing) > 1) {
    pheno_matched[non_missing] <- scale(pheno_matched[non_missing])
  }
  
  # Perform linear regression GWAS
  linear_pheno <- matrix(pheno_matched, ncol = 1)
  rownames(linear_pheno) <- fam_ids
  outfile_trans <- paste0(outfile, ".trans")
  cat("Running linear regression GWAS on transformed phenotype...\n")
  linear_gwas_parallel(filename, linear_pheno,
                       nr_blocks, outfile_trans)
  
  cat("\n")
  cat(sprintf("GWAS completed. Results written to %s.pheno1\n\n", outfile_trans))
}