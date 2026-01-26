#' Estimate genetic covariance of strata
#'
#' @useDynLib StratGWAS, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#'
#' @param strata An object returned from stratify()
#' @param filename Prefix of genotype .bed file
#' @param nr_blocks Block size for reading in genotype data (default: 1000)
#' @param outfile Name of output file (should match previous step)
#' @param SumHer Indicates whether an implementation of SumHer will be used (default: SumHer = T) or LDSC (SumHer = F)
#' @param ss_list Optional: list of length K+2, containing input summary statistics specific to strata, binary trait, and stratification variable
#' @param lds Optional: data frame containing LD scores
#' @return Returns covariance matrix of the strata
#' @export
compute_gencov <- function(strata, filename, nr_blocks = 1000, outfile, SumHer = T, ss_list = NULL, lds = NULL) {

  # Check input data
  compute_gencov_checks(strata, filename, nr_blocks = 1000, outfile, SumHer = T, ss_list = NULL, lds = NULL)

  # Read in data as multivariate phenotype
  K <- strata$K
  K_tot <- K + 2

  ids <- as.character(read.table(paste0(filename, ".fam"))[, 2])
  multi <- matrix(0, length(ids), strata$K)

  # Match phenotype data for each stratum
  for(k in 1:K) {
    multi[, k] <- strata[[paste0("group", k)]][match(ids,strata[[paste0("group", k)]][, 1]), 3]
  }

  # Scale phenotypes while preserving missing values
  for (k in 1:K) {
    non_missing <- !is.na(multi[, k])
    if (sum(non_missing) > 1) {
      multi[non_missing, k] <- scale(as.numeric(multi[non_missing, k]))
    }
  }

  rownames(multi) <- ids

  # Add binary trait and stratification variable
  y_vals <- as.numeric(strata$y[match(ids, strata$y[, 1]), 3])
  Z_vals <- as.numeric(strata$Z[match(ids, strata$Z[, 1]), 3])
  
  # Scale values
  if (sum(!is.na(y_vals)) > 1) y_vals <- as.numeric(scale(y_vals))
  if (sum(!is.na(Z_vals)) > 1) Z_vals <- as.numeric(scale(Z_vals))

  # Create new phenotype file
  multi <- cbind(multi, y_vals, Z_vals)
  multi_pheno <- cbind(ids, ids, multi)

  # Write phenotype file
  write.table(multi_pheno, paste0(outfile, ".strata"), 
              quote = FALSE, row.names = FALSE, col.names = FALSE)

  # Perform a linear regression on the data: (i) subtypes; (ii) disease; (iii) stratification variable
  if (is.null(ss_list)) {
    linear_gwas(filename, multi, nr_blocks, outfile)

    # Read in linear regression results
    ss_list <- vector("list", K_tot)
    for (k in 1:K_tot) {
      ss_list[[k]] <- read.table(paste0(outfile, ".pheno", k), header = TRUE)
    }
  }

  # Compute LD scores if not provided
  if (is.null(lds)) {
    if (length(ids) > 1000) {
      set.seed(123)  # For reproducibility
      geno_set <- sort(sample(seq_len(length(ids)), size = 1000))
    } else {
      geno_set <- seq_len(length(ids))
    }

    lds <- computeLDscoresFromBED(filename, geno_set)
    write.table(lds, paste0(outfile, ".ldscores"), 
                quote = FALSE, row.names = FALSE, col.names = TRUE)
  }

  # Initialize results matrices
  hers <- rep(NA, K_tot)
  gencov <- matrix(NA, nrow = K_tot, ncol = K_tot)
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

          cat(sprintf("SNP heritability of trait %d: %.4f (SE = %.4f)\n", 
                      i, sum$h2_snp, sum$se_h2))
        } else {
          sum_cov <- sumher_cov(ss_list[[i]], ss_list[[j]], ldscores)
          gencov[i, j] <- gencov[j, i] <- sum_cov$h2_AB

          cat(sprintf("Genetic covariance between trait %d and %d: %.4f (SE = %.4f)\n", 
                      i, j, sum_cov$h2_AB, sum_cov$se_h2_AB))
        }
      }
    }
  }

  cat("\n")

  # Compute genetic covariance matrix
  gencor <- gencov
  for (k in 1:K_tot) {
    if (is.na(hers[k]) || hers[k] <= 0) {
      if (!is.na(hers[k]) && hers[k] < 0) {
        cat(sprintf("SNP heritability of trait %d is negative, so will not compute genetic correlation\n", k))
      }
      gencor[k, ] <- NA
      gencor[, k] <- NA
    } else {
      gencor[k, ] <- gencor[k, ] / sqrt(hers[k])
      gencor[, k] <- gencor[, k] / sqrt(hers[k])
    }
  }

  # Write output files
  write.table(hers, paste0(outfile, ".hers"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(gencov, paste0(outfile, ".gencov"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(gencor, paste0(outfile, ".gencor"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)

  return(gencov)
}