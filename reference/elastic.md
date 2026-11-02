# Fit elastic net model

Development version: fits a variational Bayes elastic-net polygenic
score.

## Usage

``` r
elastic(
  filename,
  pheno_mat,
  h2 = 0.3,
  alpha_param = -0.25,
  p_en = 0.1,
  F_en = 0.5,
  chunk_size = 256,
  tol = 1e-06,
  max_scans = 10
)
```

## Arguments

- filename:

  PLINK prefix (.bed/.bim/.fam) of the genotype data

- pheno_mat:

  Numeric matrix of phenotypes (one column per phenotype), with IIDs as
  rownames

- h2:

  Assumed SNP heritability (default: 0.3)

- alpha_param:

  Selection parameter for the dependency of SNP heritability on MAF
  (default: -0.25)

- p_en, F_en:

  Elastic-net prior hyperparameters: proportion of SNPs in the
  large-effect component, and the fraction of heritability they explain
  (defaults: 0.1 and 0.5)

- chunk_size:

  Number of SNPs read per block (default: 256)

- tol:

  Convergence tolerance (default: 1e-6)

- max_scans:

  Maximum number of passes through the data (default: 10)

## Value

List returned by the C++ elastic-net solver
