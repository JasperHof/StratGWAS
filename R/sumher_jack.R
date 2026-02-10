#' Jackknife for estimating SE of SNP Heritability
#'
#' @param ss Summary statistics data.frame with columns: SNP, N, Chisq
#' @param ldscores Numeric vector of LD scores
#' @param fit_intercept Logical, whether to fit intercept (default TRUE)
#' @param tol Convergence tolerance (default 1e-6)
#' @param max_iter Maximum iterations (default 100)
#' @param alpha Selection parameter for MAF
#' @param B Number of jackknife blocks (default 100)
#' @return List with heritability estimates and likelihood
#' @export
sumher_jack <- function(ss, ldscores, 
                        fit_intercept = TRUE,
                        tol = 1e-6,
                        max_iter = 100,
                        alpha = -0.25,
                        B = 100) {

  make_blocks <- function(M, B) {
    split(seq_len(M), cut(seq_len(M), B, labels = FALSE))
  }

  M <- nrow(ss)
  blocks <- make_blocks(M, B)

  theta_jack <- numeric(B)

  for (b in seq_len(B)) {
    keep <- setdiff(seq_len(M), blocks[[b]])

    res_b <- sumher(
      ss[keep, ],
      ldscores[keep],
      fit_intercept = fit_intercept,
      alpha = alpha
    )

    theta_jack[b] <- res_b$h2_snp
  }

  theta_bar <- mean(theta_jack)
  var_jack <- (B - 1) / B * sum((theta_jack - theta_bar)^2)
  se_jack <- sqrt(var_jack)

  result <- list(
    se = se_jack,
    ests = theta_jack
  )

  return(result)
}


#' Jackknife for estimating SE of covariance
#'
#' @param ss Summary statistics data.frame with columns: SNP, N, Chisq
#' @param ldscores Numeric vector of LD scores
#' @param fit_intercept Logical, whether to fit intercept (default TRUE)
#' @param tol Convergence tolerance (default 1e-6)
#' @param max_iter Maximum iterations (default 100)
#' @param alpha Selection parameter for MAF
#' @param B Number of jackknife blocks (default 100)
#' @return List with heritability estimates and likelihood
#' @export
sumher_cov_jack <- function(ss1, ss2, ldscores,
                        fit_intercept = TRUE,
                        tol = 1e-6,
                        max_iter = 100,
                        alpha = -0.25,
                        B = 100) {

  make_blocks <- function(M, B) {
    split(seq_len(M), cut(seq_len(M), B, labels = FALSE))
  }

  M <- nrow(ss1)
  blocks <- make_blocks(M, B)

  theta_jack <- numeric(B)

  for (b in seq_len(B)) {
    keep <- setdiff(seq_len(M), blocks[[b]])

    res_b <- sumher_cov(
      ss1[keep, ],
      ss2[keep, ],
      ldscores[keep],
      fit_intercept = fit_intercept,
      alpha = alpha
    )

    theta_jack[b] <- res_b$h2_AB
  }

  theta_bar <- mean(theta_jack)
  var_jack <- (B - 1) / B * sum((theta_jack - theta_bar)^2)
  se_jack <- sqrt(var_jack)

  result <- list(
    se = se_jack,
    ests = theta_jack
  )

  return(result)
}