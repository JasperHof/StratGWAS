# Jackknife for estimating SE of SNP Heritability

Jackknife for estimating SE of SNP Heritability

## Usage

``` r
sumher_jack(
  ss,
  ldscores,
  fit_intercept = TRUE,
  tol = 1e-06,
  max_iter = 100,
  alpha = -0.25,
  B = 100
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

- alpha:

  Selection parameter for MAF

- B:

  Number of jackknife blocks (default: 100)

## Value

List with heritability estimates and likelihood
