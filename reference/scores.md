# Compute polygenic scores

Compute polygenic scores

## Usage

``` r
scores(filename, snp_input, beta_input, block_size = 256)
```

## Arguments

- filename:

  PLINK prefix (.bed/.bim/.fam) of the genotype data

- snp_input:

  Character vector of SNP IDs

- beta_input:

  Numeric vector of effect sizes, in the same order as `snp_input`

- block_size:

  Number of SNPs read per block (default: 256)

## Value

Data frame with a polygenic score per individual
