# StratGWAS

StratGWAS is an R package for modelling disease heterogeneity of binary traits in genome-wide association studies (GWAS). StratGWAS can leverage secondary phenotypes (e.g. age at diagnosis, comorbodities) as continuous and/or categorical variables to stratify cases, estimate genetic covariance between strata, and construct a transformed phenotype that maximizes power to detect genetic associations.

---

## Installation

StratGWAS can be installed from GitHub:

```r
# install.packages("devtools")
devtools::install_github("JasperHof/StratGWAS")
```

Then load the package:

```r
library(StratGWAS)
```

Please note that StratGWAS requires Rtools to install, so this should be downloaded prior to installing StratGWAS.

---

## Input Data Format

All phenotype, covariate, and stratification variables must follow **PLINK-style format**:

| Column | Description                                       |
| ------ | ------------------------------------------------- |
| 1      | Family ID (FID)                                   |
| 2      | Individual ID (IID)                               |
| 3+     | Phenotype / covariates / stratification variables |

Additional requirements:

* Binary phenotypes coded as `0 = control`, `1 = case`
* Missing values coded as `NA`
* Genotypes supplied in PLINK `.bed/.bim/.fam` format

---

## Workflow Overview

StratGWAS consists of four main steps:

1. **Stratification** (`stratify`)
   Divide cases into subgroups based on stratification variables.

2. **Genetic Covariance Estimation** (`compute_gencov`)
   Estimate genetic covariances between strata using SumHer.

3. **Phenotype Transformation** (`transform`)
   Create a transformed phenotype based on the covariance structure.

4. **GWAS Analysis** (`linear`)
   Perform GWAS on the transformed phenotype.

---

## Example Workflow

```r
library(StratGWAS)

# Load example data
data(pheno)
data(strat_cont)
data(strat_cat)
data(cov)

# Genotype file (PLINK)
filename_bed <- system.file("extdata", "data.bed", package = "StratGWAS")
filename <- gsub(".bed", "", filename_bed)

outfile <- tempfile("stratgwas_output")

# Step 1: Stratify
strata <- stratify(pheno, strat_cont = strat_cont, strat_cat = strat_cat, K = 5)

# Step 2: Genetic covariance
gencov <- compute_gencov(strata, filename, nr_blocks = 1000, outfile)

# Step 3: Transform phenotype
trans <- transform(strata, gencov, outfile)

# Step 4: Run GWAS
results <- linear(trans, filename, outfile, nr_blocks = 1000, cov = cov)

head(results$results)
```

---

## Output Files

The workflow produces several files using the provided output prefix:

### From `compute_gencov()`

* `<outfile>.strata` – Multi-phenotype file for all strata
* `<outfile>.pheno1-K` – GWAS results per stratum
* `<outfile>.ldscores` – LD scores
* `<outfile>.hers` – SNP heritabilities of strata
* `<outfile>.gencov` – Genetic covariance matrix
* `<outfile>.gencor` – Genetic correlation matrix

### From `transform()`

* `<outfile>.transformation` – Transformed phenotype file

### From `linear()`

* `<outfile>.trans.pheno1` – Final GWAS results

---

## Performance Notes

* Ensure sufficient sample size per stratum.
* Quality control of genotypes is recommended before running StratGWAS.
* StratGWAS can benefit from multithreading
