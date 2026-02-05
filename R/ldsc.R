#' LD score regression for genetic correlation
#'
#' Compute genetic correlation based on two sets of summary statistics
#'
#' @param ss1 Summary stats data.frame for trait 1
#' @param ss2 Summary stats data.frame for trait 2
#' @param ldscores Numeric vector of LD scoees
#' @return Numeric vector: h2
#' @export
ldsc_cor <- function(ss1, ss2, ldscores) {

  ldscores = as.numeric(ldscores)

  if (nrow(ss1) != length(ldscores) ||
      nrow(ss2) != length(ldscores))  stop("Summary statistics and LD scores must match in length")

  if (!all(ss1$SNP == ss2$SNP))
    warning("SNP order differs between summary statistics")
  if (nrow(ss1) < 1e4)
    warning("Very small number of SNPs for LDSC")

  # Z-scores
  z1 <- ss1$Beta / ss1$SE
  z2 <- ss2$Beta / ss2$SE
  chisq1 <- z1^2
  chisq2 <- z2^2

  M <- nrow(ss1) # Number of SNPs

  # Sample sizes
  N1 <- mean(ss1$N)
  N2 <- mean(ss2$N)
  N12 <- sqrt(N1 * N2)

  # Cross-product for genetic covariance
  y <- z1 * z2
  
  # Weights to guard against heteroscedasticity
  w <- 1 / pmax(ldscores, 1)

  # Cross-trait LDSC regression: E[z1*z2] = intercept + (N12 * cov_g / M) * ℓ
  fit_cross <- lm(y ~ ldscores, weights = w)
  beta_cross <- coef(fit_cross)["ldscores"]
  intercept_cross <- coef(fit_cross)["(Intercept)"]

  # Genetic covariance (correcting for intercept)
  cov_g <- (beta_cross * M - intercept_cross * N12) / N12

  # Heritability regressions: E[χ²] = 1 + Na + (N * h² / M) * ℓ
  fit_h2_1 <- lm(chisq1 ~ ldscores, weights = w)
  fit_h2_2 <- lm(chisq2 ~ ldscores, weights = w)
  
  beta_h2_1 <- coef(fit_h2_1)["ldscores"]
  beta_h2_2 <- coef(fit_h2_2)["ldscores"]
  intercept_h2_1 <- coef(fit_h2_1)["(Intercept)"]
  intercept_h2_2 <- coef(fit_h2_2)["(Intercept)"]

  # Heritabilities (correcting for intercept)
  h2_1 <- (beta_h2_1 * M - (intercept_h2_1 - 1) * N1) / N1
  h2_2 <- (beta_h2_2 * M - (intercept_h2_2 - 1) * N2) / N2

  # Constrain heritabilities to be non-negative
  h2_1 <- max(0, h2_1)
  h2_2 <- max(0, h2_2)

  # Genetic correlation - can only compute when h2_1 and h2_2 are positive
  rg <- NA_real_
  if (h2_1 > 0 && h2_2 > 0) {
    rg <- cov_g / sqrt(h2_1 * h2_2)
    # Constrain to valid correlation range
    rg <- max(-1, min(1, rg))
  }

  return(list(
    rg = rg,
    cov_g = cov_g,
    h2_1 = h2_1,
    h2_2 = h2_2,
    intercepts = list(
      cross = intercept_cross,
      h2_1 = intercept_h2_1,
      h2_2 = intercept_h2_2
    )
  ))
}