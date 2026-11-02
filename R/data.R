#' Example binary phenotype
#'
#' Simulated case-control phenotype for the 5,000 individuals in the example
#' genotype data (`inst/extdata/data.bed`).
#' @format A data frame with 5000 rows: `FID`, `IID`, `Phen1` (0 = control, 1 = case).
"pheno"

#' Example continuous stratification variable
#'
#' Simulated continuous secondary phenotype (e.g. age at onset) for the example data.
#' @format A data frame with 5000 rows: `FID`, `IID`, `Phen1` (numeric).
"strat_cont"

#' Example categorical stratification variable
#'
#' Simulated categorical secondary phenotype with three categories (0, 1, 2).
#' @format A data frame with 5000 rows: `FID`, `IID`, `Phen1` (integer category).
"strat_cat"

#' Example covariates
#'
#' Four simulated principal components for the example data.
#' @format A data frame with 5000 rows: `FID`, `IID`, `pc1`–`pc4`.
"cov"