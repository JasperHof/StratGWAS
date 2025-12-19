#' LD score regression for genetic correlation
#'
#' Compute genetic correlation based on two sets of summary statistics
#'
#' @param sumstats1_file Summary statistics of phenotype 1
#' @param sumstats2_file Summary statistics of phenotype 2
#' @param ldscores Numeric vector of LD scoees
#' @return Numeric vector: h2
#' @export
ldsc_cor <- function(sumstats1_file, sumstats2_file, ldscores) {
  
  # Read in summary statistics
  ss1 <- read.table(sumstats1_file, header = TRUE)
  ss2 <- read.table(sumstats2_file, header = TRUE)

  if (nrow(ss1) != length(ldscores) ||
      nrow(ss2) != length(ldscores))  stop("Summary statistics and LD scores must match in length")

  if (!all(ss1$SNP == ss2$SNP))
    warning("SNP order differs between summary statistics")

  # Z-scores
  z1 <- ss1$Beta / ss1$SE
  z2 <- ss2$Beta / ss2$SE
  y <- z1 * z2

  M <- nrow(ss1)

  # Sample sizes
  N1 <- mean(ss1$N)
  N2 <- mean(ss2$N)
  N12 <- min(N1, N2)   # conservative default

  # Cross-trait LDSC regression
  fit <- lm(y ~ ldscores)
  beta <- coef(fit)["ldscores"]

  # Genetic covariance
  cov_g <- beta * M / N12

  # --- Heritabilities (needed for rg) ---
  h2_1 <- coef(lm(ss1$Chisq - 1 ~ ldscores))["ldscores"] * M / N1
  h2_2 <- coef(lm(ss2$Chisq - 1 ~ ldscores))["ldscores"] * M / N2

  # Genetic correlation
  rg <- cov_g / sqrt(h2_1 * h2_2)

  return(list(
    rg = rg,
    cov_g = cov_g,
    h2_1 = h2_1,
    h2_2 = h2_2
  ))

  return(h2)
}