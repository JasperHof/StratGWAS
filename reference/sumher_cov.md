# Implementation of SumHer for estimating genetic covariance

Implementation of SumHer for estimating genetic covariance

## Usage

``` r
sumher_cov(
  ss1,
  ss2,
  ldscores,
  h2_1 = NULL,
  h2_2 = NULL,
  sample_overlap = 0,
  fit_intercept = TRUE,
  tol = 1e-06,
  max_iter = 100,
  alpha = -0.25
)
```

## Arguments

- ss1:

  Summary statistics for trait 1: data.frame with columns Predictor, N,
  BETA, SE, Chisq and MAF (same SNPs and order as `ldscores`)

- ss2:

  Summary statistics for trait 2, in the same format and SNP order as
  `ss1`

- ldscores:

  Numeric vector of LD scores

- h2_1:

  Heritability estimate for trait 1 (optional, will estimate if NULL)

- h2_2:

  Heritability estimate for trait 2 (optional, will estimate if NULL)

- sample_overlap:

  Proportion of sample overlap (default 0)

- fit_intercept:

  Logical, fit intercept for overlap (default TRUE)

- tol:

  Convergence tolerance (default 1e-6)

- max_iter:

  Maximum iterations (default 100)

- alpha:

  Selection parameter for the dependency of SNP heritability on MAF
  (default: -0.25)

## Value

List with genetic correlation and covariance estimates
