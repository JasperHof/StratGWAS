# Linear regression GWAS on transformed phenotype

Perform genome-wide association analysis using the transformed
phenotype. Covariates can be included and will be regressed out before
analysis. This is the fourth and final step in the StratGWAS workflow.

## Usage

``` r
linear(trans, filename, outfile, nr_blocks = 1000, cov = NULL, alpha = -0.25)
```

## Arguments

- trans:

  Object returned from [`transform`](transform.md), containing the
  transformed phenotype and associated weights

- filename:

  Prefix of genotype .bed file (without the .bed extension)

- outfile:

  Name/path prefix for output files

- nr_blocks:

  Block size for reading in genotype data (default: 1000)

- cov:

  Optional: data frame with covariates in PLINK format (FID, IID, and
  covariate columns). Missing covariate values will be imputed to their
  column means. Covariates will be regressed out of the transformed
  phenotype before GWAS analysis.

## Value

GWAS results are written to .assoc

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage with covariates
data(pheno)
data(strat_cont)
data(cov)
filename_bed <- system.file("extdata", "data.bed", package = "StratGWAS")
filename <- gsub(".bed", "", filename_bed)
outfile <- tempfile("gwas")

strata <- stratify(pheno, strat_cont = strat_cont, K = 5)
gencov <- compute_gencov(strata, filename, nr_blocks = 1000, outfile)
trans <- transform(strata, gencov)
linear(trans, filename, outfile, nr_blocks = 1000, cov = cov)

gwas_results <- read.table(paste0(outfile,".assoc"))

# View top hits
top_hits <- gwas_results$gwas_results[order(gwas_results$gwas_results$Pvalue), ]
head(top_hits, 20)
} # }
```
