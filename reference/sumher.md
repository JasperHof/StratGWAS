# Implementation of SumHer for estimating SNP Heritability

Implementation of SumHer for estimating SNP Heritability

## Usage

``` r
sumher(
  ss,
  ldscores,
  fit_intercept = TRUE,
  tol = 1e-06,
  max_iter = 100,
  alpha = -0.25
)
```

## Arguments

- ss:

  Summary statistics data.frame with columns N, Chisq and MAF (one row
  per SNP, in the same order as `ldscores`)

- ldscores:

  Numeric vector of LD scores

- fit_intercept:

  Logical, whether to fit intercept (default TRUE)

- tol:

  Convergence tolerance (default 1e-6)

- max_iter:

  Maximum iterations (default 100)

- alpha:

  Selection parameter for the dependency of SNP heritability on MAF
  (default: -0.25)

## Value

List with heritability estimates and likelihood
