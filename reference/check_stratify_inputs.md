# Check inputs for stratify function

Validates all input parameters for the stratify function

## Usage

``` r
check_stratify_inputs(pheno, strat_cont, strat_cat, K)
```

## Arguments

- pheno:

  Binary input phenotype in PLINK format

- strat_cont:

  Continuous stratification variables in PLINK format

- strat_cat:

  Categorical stratification variables in PLINK format

- K:

  Number of groups for continuous variables

- cov:

  Covariates in PLINK format

## Value

NULL if all checks pass, stops with error message otherwise
