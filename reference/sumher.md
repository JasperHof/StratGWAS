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

  Summary statistics data.frame with columns: SNP, N, Chisq

- ldscores:

  Numeric vector of LD scores

- fit_intercept:

  Logical, whether to fit intercept (default TRUE)

- tol:

  Convergence tolerance (default 1e-6)

- max_iter:

  Maximum iterations (default 100)

- weights:

  Optional weights (default: 1/ldscores)

## Value

List with heritability estimates and likelihood
