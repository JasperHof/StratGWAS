#' Estimate genetic covariance of strata
#'
#' Compute genetic covariances and correlations between case strata using 
#' SumHer. This is the second step in the StratGWAS workflow.
#' 
#' @useDynLib StratGWAS, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#'
#' @param strata An object returned from \code{\link{stratify}}
#' @param filename Prefix of genotype .bed file (without the .bed extension)
#' @param nr_blocks Block size for reading in genotype data (default: 1000)
#' @param outfile Name/path prefix for output files
#' @param lds Optional: data frame containing LD scores. The first column should contain
#'   SNP IDs, while the second column contains LD scores. If NULL, LD scores will be 
#'   computed from genotype data (using a random subset of 1000 individuals if N > 1000).
#' @param ss_list Optional: list of previously computed summary statistics for strata 
#'   in strata$multi
#' @param alpha Selection parameter for dependency MAF on SNP heritability (default: -0.25)
#' @param jack Logical specifying whether jack-knife estimates are used to compute h2 standard errors (default: FALSE)
#' @param B Number of jackknife blocks for standard error of genetic covariance (default 100)
#'
#' @return Returns a list containing:
#'   \item{gencov}{Genetic covariance matrix between strata}
#'   \item{gencor}{Genetic correlation matrix between strata}
#'   \item{hers}{Vector of SNP heritabilities for each stratum}
#'   \item{var_names}{Names of variables/strata}
#'   \item{strat_details}{Detailed stratification information from input}
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' data(pheno)
#' data(strat_cont)
#' filename_bed <- system.file("extdata", "data.bed", package = "StratGWAS")
#' filename <- gsub(".bed", "", filename_bed)
#' outfile <- tempfile("gencov")
#' 
#' strata <- stratify(pheno, strat_cont = strat_cont, K = 5)
#' gencov <- compute_gencov(strata, filename, nr_blocks = 1000, outfile)
#' 
#' # Examine results
#' print(gencov$gencov)  # Genetic covariance matrix
#' print(gencov$gencor)  # Genetic correlation matrix
#' print(gencov$hers)    # SNP heritabilities
#' 
#' 
#' }
#' 
#' @export
compute_gencov <- function(strata, filename, nr_blocks = 1000, outfile,
                           lds = NULL, ss_list = NULL,
                           alpha = -0.25, jack = FALSE, B = 50) {

  # Check input data
  # compute_gencov_checks(strata, filename, nr_blocks, outfile, SumHer, lds)

  # Get the multi phenotype matrix from strata
  multi_pheno <- strata$multi

  # Extract IDs from first two columns
  ids <- as.character(multi_pheno[, 2])

  # Get phenotype matrix (excluding first two ID columns)
  multi <- multi_pheno[, -(1:2), drop = FALSE]
  multi <- apply(multi, 2, as.numeric)
  K_tot <- ncol(multi)

  # Read genotype IDs
  fam_ids <- as.character(read.table(paste0(filename, ".fam"))[, 2])

  # Match multi to genotype file order
  multi_matched <- matrix(NA_real_, length(fam_ids), K_tot)
  rownames(multi_matched) <- fam_ids
  colnames(multi_matched) <- colnames(multi)

  for(k in 1:K_tot) {
    multi_matched[, k] <- multi[match(fam_ids, ids), k]
  }

  # Scale phenotypes while preserving missing values
  for (k in 1:K_tot) {
    non_missing <- !is.na(multi_matched[, k])
    if (sum(non_missing) > 1) {
      multi_matched[non_missing, k] <- scale(as.numeric(multi_matched[non_missing, k]))
    }
  }

  # Write phenotype file
  multi_pheno_out <- cbind(fam_ids, fam_ids, multi_matched)
  write.table(multi_pheno_out, paste0(outfile, ".strata"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)

  # Perform linear regression on all phenotypes
  if (is.null(ss_list)) {
    cat("\n")
    linear_gwas_parallel(filename, multi_matched, nr_blocks, outfile)
    cat("\n")

    # Read in linear regression results
    ss_list <- vector("list", K_tot)
    for (k in 1:K_tot) {
      cat(sprintf("Reading summary statistics of stratum %d \n", k))
      ss_list[[k]] <- read.table(paste0(outfile, ".pheno", k), header = TRUE)
    }
  }

  cat("\n")

  # Compute LD scores if not provided
  if (is.null(lds)) {
    if (length(fam_ids) > 1000) {
      set.seed(123)  # For reproducibility
      geno_set <- sort(sample(seq_len(length(fam_ids)), size = 1000))
    } else {
      geno_set <- seq_len(length(fam_ids))
    }

    lds <- computeLDscoresFromBED(filename, geno_set)
    write.table(lds, paste0(outfile, ".ldscores"), 
                quote = FALSE, row.names = FALSE, col.names = TRUE)
  }
  cat("\n")

  # Initialize results matrices
  hers <- rep(NA, K_tot)
  hers_se <- rep(NA, K_tot)

  gencov <- matrix(NA, nrow = K_tot, ncol = K_tot)
  rownames(gencov) <- colnames(gencov) <- colnames(multi)
  ldscores <- lds$Tagging[match(ss_list[[1]]$Predictor, lds$Predictor)]

  # Use SumHer implementation for genetic correlations
  jack_ests <- vector("list", B)
  for (b in 1:B) jack_ests[[b]] = matrix(NA, K_tot, K_tot)

  for (i in 1:K_tot) {
    for (j in i:K_tot) {
      if (i == j) {
        sum <- sumher(ss_list[[i]], ldscores, alpha = alpha)

        if (jack == TRUE) {
          jack_sum <- sumher_jack(ss_list[[i]], ldscores, alpha = alpha, B = B)
          for (b in 1:B) jack_ests[[b]][i, j] <- jack_ests[[b]][j, i] <- jack_sum$ests[b]
    
          gencov[i, j] <- sum$h2_snp
          hers[i] <- sum$h2_snp
          hers_se[i] <- jack_sum$se
    
          cat(sprintf("SNP heritability of %s: %.4f (SE = %.4f)\n",
                      colnames(multi)[i], sum$h2_snp, jack_sum$se))
        } else {
          gencov[i, j] <- sum$h2_snp
          hers[i] <- sum$h2_snp
          hers_se[i] <- sum$se_h2
    
          cat(sprintf("SNP heritability of %s: %.4f (SE = %.4f)\n",
                      colnames(multi)[i], sum$h2_snp, sum$se_h2))
        }
      } else {
        sum_cov <- sumher_cov(ss_list[[i]], ss_list[[j]], ldscores, alpha = alpha)

        if (jack == TRUE) {
          jack_cov <- sumher_cov_jack(ss_list[[i]], ss_list[[j]], ldscores, alpha = alpha, B = B)
          for (b in 1:B) jack_ests[[b]][i, j] <- jack_ests[[b]][j, i] <- jack_cov$ests[b]
    
          gencov[i, j] <- gencov[j, i] <- sum_cov$h2_AB
          cat(sprintf("Genetic covariance between %s and %s: %.4f (SE = %.4f)\n", 
                    colnames(multi)[i], colnames(multi)[j], sum_cov$h2_AB, jack_cov$se))
        } else {
          gencov[i, j] <- gencov[j, i] <- sum_cov$h2_AB
          cat(sprintf("Genetic covariance between %s and %s: %.4f (SE = %.4f)\n", 
                    colnames(multi)[i], colnames(multi)[j], sum_cov$h2_AB, sum_cov$se_h2_AB))
        }
      }
    }
  }

  if (jack == TRUE) {
    for (b in 1:B){
      rownames(jack_ests[[b]]) <- colnames(jack_ests[[b]]) <- colnames(multi)
    }
  }

  cat("\n")

  # Compute genetic correlation matrix
  gencor <- gencov
  rownames(gencor) <- colnames(gencor) <- colnames(multi)
  
  for (k in 1:K_tot) {
    if (is.na(hers[k]) || hers[k] <= 0) {
      if (!is.na(hers[k]) && hers[k] < 0) {
        cat(sprintf("SNP heritability of %s is negative, so will not compute genetic correlation\n", 
                    colnames(multi)[k]))
      }
      gencor[k, ] <- NA
      gencor[, k] <- NA
    } else {
      gencor[k, ] <- gencor[k, ] / sqrt(hers[k])
      gencor[, k] <- gencor[, k] / sqrt(hers[k])
    }
  }

  # Compute liability scale heritabilities + SE
  prevs <- colMeans(multi_matched > 0, na.rm = T)
  t <- qnorm(1 - prevs)
  z <- dnorm(t)

  hers_liab <- hers * (prevs * (1 - prevs)) / z^2
  hers_liab_se <- hers_se * (prevs * (1 - prevs)) / z^2
  hers_all <- data.frame("h2_obs" = hers, "h2_obs_SE" = hers_se, "h2_liab" = hers_liab, "h2_liab_SE" = hers_liab_se)
  hers_all[rownames(hers_all) %in% paste0("continuous_", 1:100), c("h2_liab", "h2_liab_SE")] <- NA

  # Write output files with row/column names
  write.table(hers_all, paste0(outfile, ".hers"),
              quote = FALSE, row.names = TRUE, col.names = TRUE)
  write.table(gencov, paste0(outfile, ".gencov"),
              quote = FALSE, row.names = TRUE, col.names = TRUE)
  write.table(gencor, paste0(outfile, ".gencor"),
              quote = FALSE, row.names = TRUE, col.names = TRUE)

  # Return results with metadata
  result <- list(
    gencov = gencov,
    gencor = gencor,
    jack_ests = jack_ests,
    jack = jack,
    hers = hers,
    var_names = colnames(multi),
    strat_details = strata$strat_details
  )

  return(result)
}