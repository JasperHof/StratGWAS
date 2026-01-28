#' Estimate genetic covariance of multiple strata
#'
#' @useDynLib StratGWAS, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#'
#' @param strata An object returned from stratify_multi()
#' @param filename Prefix of genotype .bed file
#' @param nr_blocks Block size for reading in genotype data (default: 1000)
#' @param outfile Name of output file (should match previous step)
#' @param SumHer Indicates whether an implementation of SumHer will be used (default: SumHer = T) or LDSC (SumHer = F)
#' @param ss_list Optional: list of length K_tot, containing input summary statistics
#' @param lds Optional: data frame containing LD scores
#' @return Returns covariance matrix of the strata
#' @export
compute_gencov_multi <- function(strata, filename, nr_blocks = 1000, outfile, 
                                 SumHer = T, ss_list = NULL, lds = NULL) {

  # Check input data
  # compute_gencov_checks(strata, filename, nr_blocks, outfile, SumHer, ss_list, lds)

  # Get the multi phenotype matrix from strata
  multi_pheno <- strata$multi
  
  # Extract IDs from first two columns
  ids <- as.character(multi_pheno[, 1])
  
  # Get phenotype matrix (excluding first two ID columns)
  multi <- multi_pheno[, -(1:2), drop = FALSE]
  K_tot <- ncol(multi)
  
  # Read genotype IDs
  fam_ids <- as.character(read.table(paste0(filename, ".fam"))[, 2])
  
  # Match multi to genotype file order
  multi_matched <- matrix(NA, length(fam_ids), K_tot)
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
    linear_gwas(filename, multi_matched, nr_blocks, outfile)

    # Read in linear regression results
    ss_list <- vector("list", K_tot)
    for (k in 1:K_tot) {
      ss_list[[k]] <- read.table(paste0(outfile, ".pheno", k), header = TRUE)
    }
  }

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

  # Initialize results matrices
  hers <- rep(NA, K_tot)
  gencov <- matrix(NA, nrow = K_tot, ncol = K_tot)
  rownames(gencov) <- colnames(gencov) <- colnames(multi)
  ldscores <- lds$Tagging[match(ss_list[[1]]$Predictor, lds$Predictor)]

  if(SumHer == F){ 
    # Use LDSC implementation for genetic correlations
    for (i in 1:K_tot) {
      for (j in i:K_tot) {
        if (i == j) {
          ldsc <- ldsc_cor(ss_list[[i]], ss_list[[j]], ldscores)
          gencov[i, j] <- ldsc$cov_g
          hers[i] <- ldsc$h2_1
        } else {
          ldsc <- ldsc_cor(ss_list[[i]], ss_list[[j]], ldscores)
          gencov[i, j] <- gencov[j, i] <- ldsc$cov_g
        }
      }
    }
  } else { 
    # Use SumHer implementation for genetic correlations (default)
    for (i in 1:K_tot) {
      for (j in i:K_tot) {
        if (i == j) {
          sum <- sumher(ss_list[[i]], ldscores)
          gencov[i, j] <- sum$h2_snp
          hers[i] <- sum$h2_snp

          cat(sprintf("SNP heritability of %s: %.4f (SE = %.4f)\n", 
                      colnames(multi)[i], sum$h2_snp, sum$se_h2))
        } else {
          sum_cov <- sumher_cov(ss_list[[i]], ss_list[[j]], ldscores)
          gencov[i, j] <- gencov[j, i] <- sum_cov$h2_AB

          cat(sprintf("Genetic covariance between %s and %s: %.4f (SE = %.4f)\n", 
                      colnames(multi)[i], colnames(multi)[j], sum_cov$h2_AB, sum_cov$se_h2_AB))
        }
      }
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

  # Write output files with row/column names
  write.table(cbind(rownames = colnames(multi), hers), paste0(outfile, ".hers"),
              quote = FALSE, row.names = FALSE, col.names = TRUE)
  write.table(gencov, paste0(outfile, ".gencov"),
              quote = FALSE, row.names = TRUE, col.names = TRUE)
  write.table(gencor, paste0(outfile, ".gencor"),
              quote = FALSE, row.names = TRUE, col.names = TRUE)

  # Return results with metadata
  result <- list(
    gencov = gencov,
    gencor = gencor,
    hers = hers,
    var_names = colnames(multi),
    strat_details = strata$strat_details
  )
  
  return(result)
}