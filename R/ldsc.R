#' LD score regression
#'
#' Compute SNP heritability based on summary statistics
#'
#' @param sumstats_file Summary statistics
#' @param ldscores Numeric vector of LD scoees
#' @return Numeric vector: h2
#' @export
ldsc <- function(sumstats_file, ldscores) {
  
  # Read in summary statistics
  ss <- read.table(sumstats_file, head = T)

  if(dim(ss)[1] != length(ldscores)) stop("summary statistics and LD scores should be of the same size")

  y <- ss$Chisq - 1
  M <- nrow(ss)
  N <- mean(ss$N)

  # fit the regression
  fit <- lm(y ~ ldscores)
  beta <- coef(fit)["ldscores"]
  h2 <- beta * M / N

  return(h2)
}