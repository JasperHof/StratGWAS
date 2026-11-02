# Jackknife for estimating SE of covariance

Jackknife for estimating SE of covariance

## Usage

``` r
sumher_cov_jack(
  ss1,
  ss2,
  ldscores,
  fit_intercept = TRUE,
  tol = 1e-06,
  max_iter = 100,
  alpha = -0.25,
  B = 100
)
```

## Arguments

- ldscores:

  Numeric vector of LD scores

- fit_intercept:

  Logical, whether to fit intercept (default TRUE)

- tol:

  Convergence tolerance (default 1e-6)

- max_iter:

  Maximum iterations (default 100)

- alpha:

  Selection parameter for MAF

- B:

  Number of jackknife blocks (default 100)

- ss:

  Summary statistics data.frame with columns: SNP, N, Chisq

## Value

List with heritability estimates and likelihood
