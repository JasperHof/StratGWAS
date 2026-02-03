#' StratGWAS: Stratified Genome-Wide Association Studies
#'
#' @description
#' StratGWAS performs stratified genome-wide association studies (GWAS) by
#' leveraging genetic heterogeneity within cases. The package stratifies cases
#' based on continuous and/or categorical variables (e.g., age at diagnosis,
#' disease subtype), estimates genetic covariances between strata, and creates
#' a transformed phenotype that maximizes power to detect genetic associations.
#'
#' @details
#' The StratGWAS workflow consists of four main steps:
#'
#' \strong{1. Stratification} (\code{\link{stratify}}):
#' Divide cases into subgroups based on stratification variables
#'
#' \strong{2. Genetic Covariance Estimation} (\code{\link{compute_gencov}}):
#' Estimate genetic covariances between strata using SumHer or LDSC
#'
#' \strong{3. Phenotype Transformation} (\code{\link{transform}}):
#' Create a transformed phenotype based on the genetic covariance structure
#'
#' \strong{4. GWAS Analysis} (\code{\link{linear}}):
#' Perform genome-wide association analysis on the transformed phenotype
#'
#' @section Key Functions:
#' \describe{
#'   \item{\code{\link{stratify}}}{Create case strata from continuous/categorical variables}
#'   \item{\code{\link{compute_gencov}}}{Estimate genetic covariance matrix between strata}
#'   \item{\code{\link{transform}}}{Generate transformed phenotype for GWAS}
#'   \item{\code{\link{linear}}}{Run linear regression GWAS on transformed phenotype}
#' }
#'
#' @section Input Data Format:
#' All phenotype, covariate and stratification files should be in PLINK format:
#' \itemize{
#'   \item Column 1: Family ID (FID)
#'   \item Column 2: Individual ID (IID)
#'   \item Column 3+: Phenotype or stratification variables
#' }
#'
#' Binary phenotypes should be coded as 0 (controls) and 1 (cases).
#' Missing values should be coded as NA.
#'
#' @examples
#' \dontrun{
#' # ============================================================================
#' # EXAMPLE 1: Basic workflow with continuous stratification variable
#' # ============================================================================
#'
#' # Load example data
#' data(pheno)           # Binary phenotype (case-control)
#' data(strat_cont)      # Continuous stratification variable (e.g., age at diagnosis)
#' data(strat_cat)       # Categorical stratification variable (e.g., comorbidity)
#' data(cov)             # Covariates (optional)
#'
#' # Specify genotype file (PLINK .bed format)
#' filename_bed <- system.file("extdata", "data.bed", package = "StratGWAS")
#' filename <- gsub(".bed", "", filename_bed)
#'
#' # Specify output file prefix
#' outfile <- tempfile("stratgwas_output")
#'
#' # Step 1: Stratify cases into 5 quantile groups
#' strata <- stratify(pheno, strat_cont = strat_cont, strat_cat = strat_cat, K = 5)
#'
#' # Examine stratification results
#' table(strata$info$groups)           # Number of cases per stratum
#' length(strata$strat_miss)           # Cases with missing stratification values
#' head(strata$multi)                  # Multi-phenotype matrix
#' 
#' # Step 2: Compute genetic covariances between strata
#' gencov <- compute_gencov(strata, filename, nr_blocks = 1000, outfile)
#' 
#' # Examine genetic covariance results
#' print(gencov$gencov)                # Genetic covariance matrix
#' print(gencov$gencor)                # Genetic correlation matrix
#' print(gencov$hers)                  # SNP heritabilities
#' 
#' # Step 3: Transform phenotype based on genetic covariance structure
#' trans <- transform(strata, gencov)
#' 
#' # Examine transformation results
#' head(trans$transformed_pheno)       # Transformed phenotype values
#' print(trans$weights)                # Transformation weights
#' print(trans$inflation_criteria)     # Expected inflation for continuous variables
#' 
#' # Step 4: Run GWAS on transformed phenotype with covariates
#' results <- linear(trans, filename, outfile, nr_blocks = 1000, cov = cov)
#' 
#' # Examine GWAS results
#' head(results$results)               # Top SNP associations
#' results$n_samples                   # Sample size used
#' 
#' # ============================================================================
#' # EXAMPLE 2: Categorical stratification (e.g., disease subtypes)
#' # ============================================================================
#' 
#' data(strat_cat)       # Categorical variable with multiple subtypes
#' 
#' # Stratify by categorical variable
#' strata_cat <- stratify(pheno, strat_cat = strat_cat)
#' 
#' # Check detected categories
#' strata_cat$strat_details[[1]]$categories
#' strata_cat$strat_details[[1]]$K
#' 
#' # Continue with genetic covariance estimation
#' gencov_cat <- compute_gencov(strata_cat, filename, nr_blocks = 1000, outfile)
#' trans_cat <- transform(strata_cat, gencov_cat)
#' results_cat <- linear(trans_cat, filename, outfile, nr_blocks = 1000, cov = cov)
#' 
#' # ============================================================================
#' # EXAMPLE 3: Multiple stratification variables
#' # ============================================================================
#' 
#' # Stratify by both continuous and categorical variables
#' strata_multi <- stratify(pheno,
#'                          strat_cont = strat_cont,   # e.g., age at diagnosis
#'                          strat_cat = strat_cat,     # e.g., disease subtype
#'                          K = 5)
#' 
#' # Check stratification structure
#' strata_multi$n_cont                 # Number of continuous variables
#' strata_multi$n_bin                  # Number of categorical variables
#' ncol(strata_multi$multi) - 3        # Total number of strata + variables
#' 
#' # Complete workflow
#' gencov_multi <- compute_gencov(strata_multi, filename, nr_blocks = 1000, outfile)
#' trans_multi <- transform(strata_multi, gencov_multi)
#' results_multi <- linear(trans_multi, filename, outfile, nr_blocks = 1000, cov = cov)
#' 
#' # ============================================================================
#' # EXAMPLE 4: Adjust number of strata for continuous variables
#' # ============================================================================
#' 
#' # Use 10 quantile groups instead of 5
#' strata_10 <- stratify(pheno, strat_cont = strat_cont, K = 10)
#' gencov_10 <- compute_gencov(strata_10, filename, nr_blocks = 1000, outfile)
#' trans_10 <- transform(strata_10, gencov_10)
#' results_10 <- linear(trans_10, filename, outfile, nr_blocks = 1000, cov = cov)
#' 
#' # ============================================================================
#' # EXAMPLE 5: GWAS without covariates
#' # ============================================================================
#' 
#' results_nocov <- linear(trans, filename, outfile, nr_blocks = 1000)
#' 
#' # ============================================================================
#' # EXAMPLE 6: Using LDSC instead of SumHer
#' # ============================================================================
#' 
#' gencov_ldsc <- compute_gencov(strata, filename, nr_blocks = 1000,
#'                               outfile, SumHer = FALSE)
#' trans_ldsc <- transform(strata, gencov_ldsc)
#' results_ldsc <- linear(trans_ldsc, filename, outfile, nr_blocks = 1000, cov = cov)
#' 
#' # ============================================================================
#' # OUTPUT FILES
#' # ============================================================================
#' 
#' # The workflow generates several output files with the specified prefix:
#' # 
#' # From compute_gencov():
#' #   <outfile>.strata      - Multi-phenotype file for all strata
#' #   <outfile>.pheno1-K    - GWAS results for each stratum
#' #   <outfile>.ldscores    - LD scores
#' #   <outfile>.hers        - SNP heritabilities
#' #   <outfile>.gencov      - Genetic covariance matrix
#' #   <outfile>.gencor      - Genetic correlation matrix
#' #
#' # From linear():
#' #   <outfile>.transformed_pheno - Transformed phenotype file
#' #   <outfile>.pheno1            - Final GWAS results
#' }
#'
#' @author Jasper Hof \email{jasper.hof@qgg.au.dk}
#'
#' @keywords internal
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
#' @useDynLib StratGWAS, .registration = TRUE
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL