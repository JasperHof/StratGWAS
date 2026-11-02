#' @useDynLib StratGWAS, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom RcppParallel setThreadOptions
#' @importFrom stats dnorm lm pchisq predict qnorm quantile residuals rnorm smooth.spline
#' @importFrom utils read.table write.table
NULL