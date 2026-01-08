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

  M <- nrow(ss1)

  # Sample sizes
  N1 <- mean(ss1$N)
  N2 <- mean(ss2$N)
  N12 <- sqrt(N1 * N2)

  y <- z1 * z2
  w <- 1 / pmax(ldscores, 1)    # Weights to guard against heteroscadasticy

  # Cross-trait LDSC regression
  fit <- lm(y ~ ldscores, weights = w)
  beta <- coef(fit)["ldscores"]
  intercept <- coef(fit)[1]

  # Genetic covariance
  cov_g <- beta * M / N12

  # Heritabilities
  ss1$Chisq_1 = ss1$Chisq - 1
  ss2$Chisq_1 = ss2$Chisq - 1

  fit_h2_1 <- lm(ss1$Chisq ~ ldscores, weights = w)
  fit_h2_2 <- lm(ss2$Chisq ~ ldscores, weights = w)
  h2_1 <- coef(fit_h2_1)["ldscores"] * M / N1
  h2_2 <- coef(fit_h2_2)["ldscores"] * M / N2
  intercept_h2_1 <- coef(fit_h2_1)[1]
  intercept_h2_2 <- coef(fit_h2_2)[1]

  # Genetic correlation - can only compute when h1 and h2 are positive
  rg <- cov_g / sqrt(max(h2_1, 0.001) * max(h2_2, 0.001))

  rg <- NA_real_
  if (!is.na(h2_1) && !is.na(h2_2) && h2_1 > 0 && h2_2 > 0) {
    rg <- cov_g / sqrt(h2_1 * h2_2)
  }

  return(list(
    rg = rg,
    cov_g = cov_g,
    h2_1 = h2_1,
    h2_2 = h2_2,
    intercepts = list(
      cross = intercept,
      h2_1 = intercept_h2_1,
      h2_2 = intercept_h2_2
    )
  ))
}