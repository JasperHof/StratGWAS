# Estimate genetic covariance of strata

Compute genetic covariances and correlations between case strata using
SumHer. This is the second step in the StratGWAS workflow.

## Usage

``` r
compute_gencov(
  strata,
  filename,
  nr_blocks = 1000,
  outfile,
  cov = NULL,
  lds = NULL,
  ss_list = NULL,
  alpha = -0.25,
  jack = TRUE,
  B = 50
)
```

## Arguments

- strata:

  An object returned from [`stratify`](stratify.md)

- filename:

  Prefix of genotype .bed file (without the .bed extension)

- nr_blocks:

  Block size for reading in genotype data (default: 1000)

- outfile:

  Name/path prefix for output files

- cov:

  Optional: covariates in PLINK format (FID, IID, then covariates).
  Missing values are imputed by the column mean, and covariates are
  regressed out of each stratum before the GWAS.

- lds:

  Optional: data frame containing LD scores. The first column should
  contain SNP IDs, while the second column contains LD scores. If NULL,
  LD scores will be computed from genotype data (using a random subset
  of 1000 individuals if N \> 1000).

- ss_list:

  Optional: list of previously computed summary statistics for strata in
  strata\$multi

- alpha:

  Selection parameter for dependency MAF on SNP heritability (default:
  -0.25)

- jack:

  Logical specifying whether jack-knife estimates are used to compute h2
  standard errors (default: FALSE)

- B:

  Number of jackknife blocks for standard error of genetic covariance
  (default: 50)

## Value

Returns a list containing:

- gencov:

  Genetic covariance matrix between strata

- gencor:

  Genetic correlation matrix between strata

- hers:

  Vector of SNP heritabilities for each stratum

- var_names:

  Names of variables/strata

- strat_details:

  Detailed stratification information from input

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage
data(pheno)
data(strat_cont)
filename_bed <- system.file("extdata", "data.bed", package = "StratGWAS")
filename <- gsub(".bed", "", filename_bed)
outfile <- tempfile("gencov")

strata <- stratify(pheno, strat_cont = strat_cont, K = 5)
gencov <- compute_gencov(strata, filename, nr_blocks = 1000, outfile)

# Examine results
print(gencov$gencov)  # Genetic covariance matrix
print(gencov$gencor)  # Genetic correlation matrix
print(gencov$hers)    # SNP heritabilities


} # }
```
