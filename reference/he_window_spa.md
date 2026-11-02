# SPA solver for HE-based variance components with correction

Development version of [`sparrow`](sparrow.md); prefer
[`sparrow()`](sparrow.md) for new analyses.

## Usage

``` r
he_window_spa(
  filename,
  pheno,
  out_file,
  annotation = NULL,
  annot_names = NULL,
  loco_prs = NULL,
  prs_mask = NULL,
  chunk_size = 256,
  window_bp = 1e+06,
  do_windows = T,
  alpha = -1,
  covariates = NULL,
  cov_df = NA,
  SPA = T,
  spa_pval_threshold = 0.1,
  binary = F,
  binary_raw = NULL,
  coher = F,
  chr = NULL,
  batch_size = 64,
  n_threads = 0,
  max_window_chunks = 8,
  genes = NULL,
  kvik_summaries = NULL,
  acatv_mac = 10
)
```

## Arguments

- filename:

  PLINK prefix (.bed/.bim/.fam) holding the variants

- pheno:

  Data frame: FID, IID, then one column per phenotype

- out_file:

  Output prefix

- annotation:

  Optional 0/1 integer matrix, one row per .bim variant and one column
  per category

- annot_names:

  Optional category names

- loco_prs:

  Optional matrix of LDAK-KVIK Step 1 `.loco.prs` files

- prs_mask:

  Optional categories x PRS sources matrix: which full PRS is subtracted

- chunk_size, window_bp, max_window_chunks:

  Region and window sizes

- do_windows:

  Also estimate ~1 Mb windows (default: TRUE)

- alpha:

  Heritability model(s) tested (default: -1)

- covariates:

  Optional numeric covariate matrix

- cov_df:

  Degrees of freedom used by covariates (NA = determined automatically)

- SPA:

  Compute saddlepoint p-values (default: TRUE)

- spa_pval_threshold:

  Apply SPA only to tests with p-values below this threshold

- binary:

  TRUE if the phenotypes are binary

- binary_raw:

  Optional matrix of the raw 0/1 phenotypes (for binary traits)

- coher:

  TRUE to estimate co-heritability for all pairs of phenotypes

- chr:

  Optional single chromosome to analyse

- batch_size:

  Units per parallel batch

- n_threads:

  Number of threads (0 = RcppParallel default)

- genes:

  Optional gene definitions: gene, chr, start, end

- kvik_summaries:

  Optional LDAK-KVIK Step 2 `.summaries` files

- acatv_mac:

  Variants with MAC below this are collapsed into a burden in ACAT-V

## Value

NULL; results are written to `out_file`
