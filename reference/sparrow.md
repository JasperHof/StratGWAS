# SPARROW: region- and gene-based (co-)heritability analysis of rare variants

Estimates the heritability explained by genes / regions (per variant
category) and tests it using a saddlepoint approximation. Phenotypes are
matched to the genotypes by IID, covariates are regressed out of the
phenotypes here (and out of the genotypes inside the C++ code), and all
bookkeeping that used to be manual (cov_df, binary_raw, the logistic
null model, prs_mask) is done automatically.

## Usage

``` r
sparrow(
  filename,
  pheno,
  out_file,
  covariates = NULL,
  annotation = NULL,
  annot_names = NULL,
  genes = NULL,
  alpha = c(-1, -0.5, 0),
  loco_prs = NULL,
  kvik_summaries = NULL,
  binary = FALSE,
  coher = FALSE,
  h2_without_prs = TRUE,
  prs_mask = NULL,
  do_windows = TRUE,
  chunk_size = 256,
  window_bp = 1e+06,
  max_window_chunks = 8,
  max_gene_snps = 512,
  min_gene_snps = 2,
  acatv_mac = 10,
  collapse_mac = 0,
  spa = TRUE,
  chr = NULL,
  n_threads = 0,
  batch_size = 16
)
```

## Arguments

- filename:

  PLINK prefix (.bed/.bim/.fam) holding rare variants and common SNPs.

- pheno:

  Data frame or file: FID, IID, then one column per phenotype.

- out_file:

  Output prefix. Writes \<out_file\> (regions), \<out_file\>.genes and
  \<out_file\>.windows.

- covariates:

  Optional data frame or file: FID, IID, then numeric covariates (an
  intercept is added automatically; do not include one).

- annotation:

  Optional 0/1 matrix or file (with header), one column per variant
  category (e.g. rare / common) and one row per .bim variant. A file may
  start with a column of variant IDs, which is matched to the .bim. NULL
  = one category with all variants.

- annot_names:

  Category names (default: column names of `annotation`).

- genes:

  Optional data frame or file: gene, chr, start, end (extra columns
  ignored). Must be on the same genome build as the .bim.

- alpha:

  Heritability models tested; `alpha[1]` is used for the estimates.

- loco_prs:

  Optional LDAK-KVIK Step 1 `.loco.prs` files: one per phenotype
  (character vector), or a matrix with one row per phenotype and one
  column per source.

- kvik_summaries:

  Optional LDAK-KVIK Step 2 `.summaries` files, one per phenotype.

- binary:

  TRUE if the phenotypes are 0/1 (1/2 coding is converted).

- coher:

  TRUE to estimate co-heritability for all pairs of phenotypes.

- h2_without_prs:

  With `loco_prs`: estimate h2 from y (TRUE) and use y - LOCO only for
  p-values. Recommended; see the supplement.

- prs_mask:

  Advanced. Categories x PRS sources; which full PRS is subtracted in
  the window pass. Default: 0 if h2_without_prs, else 1.

- do_windows:

  Estimate ~1 Mb windows and report window-calibrated h2 (h2_adj).

- chunk_size, window_bp, max_window_chunks:

  Region and window sizes.

- max_gene_snps, min_gene_snps:

  Genes are split into blocks of at most `max_gene_snps` variants; genes
  with fewer than `min_gene_snps` are skipped.

- acatv_mac:

  Variants with MAC below this are collapsed into a burden in ACAT-V.

- collapse_mac:

  Variance-component test: variants with MAC \<= collapse_mac are summed
  into one pseudo-variant per category of each gene (as SAIGE-GENE+). 0
  = off (default). ACAT-V still uses all variants individually.

- spa:

  Compute saddlepoint p-values (FALSE = estimates only, much faster).

- chr:

  Optional single chromosome to analyse.

- n_threads:

  Number of threads (0 = RcppParallel default).

- batch_size:

  Units per parallel batch (affects progress reporting only).

## Value

Invisibly, a list with the output file names, the analysed traits and
number of individuals, and the summary returned by the C++ code.

## Examples

``` r
# \donttest{
# example data shipped with the package
dir <- system.file("extdata", package = "StratGWAS")

res <- sparrow(filename       = file.path(dir, "data"),
               pheno          = file.path(dir, "pheno_quant"),
               out_file       = file.path(tempdir(), "sparrow_example"),
               annotation     = file.path(dir, "annot"),
               genes          = file.path(dir, "genes"),
               loco_prs       = file.path(dir, "kvik.step1.loco.prs"),
               kvik_summaries = file.path(dir, "kvik.step2.summaries"),
               n_threads      = 2)
#> 5000 of 5000 individuals in `pheno` are in the .fam file
#> Dropping 0 individuals with missing phenotypes
#> SPARROW: 5000 individuals, 1 phenotype(s), 1 covariate df, alpha = -1, -0.5, 0
#> StratGWAS C++ build: Nov  2 2026 13:26:04   (p-values: one-sided VC / two-sided coher, exact binary SPA + lattice, binary_pi, ACAT-V via KVIK, h2 without PRS, collapse_mac)
#> Genotype file: 5000 individuals, 9071 SNPs
#> Fast .bed reader calibrated for genotype (codes 0/1/2/3 -> 2/-1/1/0; verified on 3 extra blocks)
#> Categories: pLoF=3098 missense=2935 synonymous=3038
#> Individuals with complete data: 5000
#> Multi-alpha: association p-values are the ACAT combination over alpha = {-1, -0.5, 0}; vg / se / h2 use alpha = -1.
#> Gene definitions: 114 read. NOTE: these coordinates must be on the same genome build as the .bim.
#> Variance-component p-values: one-sided (upper tail).
#> Environment moment: T(env,env) = n - 1
#> LOCO PRS /home/runner/work/_temp/Library/StratGWAS/extdata/kvik.step1.loco.prs: 10 chromosome columns, 5000 rows matched, 10 columns map onto the .bim chromosomes.
#>   full PRS = sum(LOCO)/9;  sd(full) per column: 0.582672
#> LOCO PRS used for p-values only: vg / se_vg / h2 of chunks and genes are estimated on y itself
#>   (an in-sample PRS shrinks heritability estimates). Use prs_mask = 0 to estimate windows on y as well.
#> PRS mask (pass 2, category x source):
#>   pLoF: 0
#>   missense: 0
#>   synonymous: 0
#> KVIK summaries /home/runner/work/_temp/Library/StratGWAS/extdata/kvik.step2.summaries: 9070 variants, 9070 matched to the .bim
#> ACAT-V: KVIK single-variant tests, weights [2f(1-f)]^(1+alpha) for each alpha, variants with MAC < 10 collapsed into a burden.
#>   p_spa = ACAT over 6 p-values (3 variance-component + 3 ACAT-V).
#> Pass 1: 40 chunks of 256 SNPs, no flanks, 3 categories, phenotype = y - LOCO PRS of the tested chromosome.
#>   pass 1: K up to ~256, ~9 MB per thread, ~0.0183258 GB at 2 threads.
#> Pass 2: 10 windows of ~1000000 bp (at most 8 chunks) with neighbouring windows as flanks, phenotype = y - full PRS per prs_mask.
#>   mean K ~ 907, worst case 908; cost per window scales as K^2.
#>   pass 2: K up to ~908, ~73 MB per thread, ~0.144395 GB at 2 threads.
#> Genes: 113 with >= 2 SNPs -> 113 blocks (split at 512 SNPs), no flanks, phenotype = y - LOCO PRS.
#>   1 genes skipped: fewer than 2 SNPs in the genotype file.
#> 
[pass1] 4/40  10%  elapsed 0s  eta 1s     
[pass1] 8/40  20%  elapsed 0s  eta 1s     
[pass1] 12/40  30%  elapsed 0s  eta 1s     
[pass1] 16/40  40%  elapsed 0s  eta 1s     
[pass1] 20/40  50%  elapsed 0s  eta 0s     
[pass1] 24/40  60%  elapsed 1s  eta 0s     
[pass1] 28/40  70%  elapsed 1s  eta 0s     
[pass1] 32/40  80%  elapsed 1s  eta 0s     
[pass1] 36/40  90%  elapsed 1s  eta 0s     
[pass1] 40/40  100%  elapsed 1s  eta 0s     
#> 
[genes] 11/113  9%  elapsed 0s  eta 0s     
[genes] 24/113  21%  elapsed 0s  eta 0s     
[genes] 35/113  30%  elapsed 0s  eta 0s     
[genes] 46/113  40%  elapsed 0s  eta 0s     
[genes] 57/113  50%  elapsed 0s  eta 0s     
[genes] 68/113  60%  elapsed 0s  eta 0s     
[genes] 81/113  71%  elapsed 0s  eta 0s     
[genes] 91/113  80%  elapsed 0s  eta 0s     
[genes] 102/113  90%  elapsed 0s  eta 0s     
[genes] 113/113  100%  elapsed 0s  eta 0s     
#> 
[pass2] 4/10  40%  elapsed 0s  eta 0s     
[pass2] 8/10  80%  elapsed 0s  eta 0s     
[pass2] 10/10  100%  elapsed 0s  eta 0s     
#> Chunks tested: 40
#> Window x category estimates: 30
#> Genes tested: 113 of 113

genes <- read.table(res$files["genes"], header = TRUE)
head(genes[genes$category == "ACAT", ])          # one p-value per gene
#>     gene chr start end block n_blocks category m_cat phenotype vg se_vg h2
#> 4  GENE1   1    25 111    NA        1     ACAT    87     Phen1 NA    NA NA
#> 8  GENE2   1   128 147    NA        1     ACAT    20     Phen1 NA    NA NA
#> 12 GENE3   1   159 220    NA        1     ACAT    62     Phen1 NA    NA NA
#> 16 GENE4   1   244 321    NA        1     ACAT    78     Phen1 NA    NA NA
#> 20 GENE5   1   350 389    NA        1     ACAT    40     Phen1 NA    NA NA
#> 24 GENE6   1   409 481    NA        1     ACAT    73     Phen1 NA    NA NA
#>    vg_env       p_spa spa_used p_alpha_.1 p_alpha_.0.5 p_alpha_0 p_vc p_acatv
#> 4      NA 3.73926e-06       NA         NA           NA        NA   NA      NA
#> 8      NA 1.75149e-05       NA         NA           NA        NA   NA      NA
#> 12     NA 1.82862e-04       NA         NA           NA        NA   NA      NA
#> 16     NA 5.07700e-02       NA         NA           NA        NA   NA      NA
#> 20     NA 6.30794e-01       NA         NA           NA        NA   NA      NA
#> 24     NA 8.13577e-01       NA         NA           NA        NA   NA      NA
#>    p_acatv_.1 p_acatv_.0.5 p_acatv_0 vg_adj h2_adj
#> 4          NA           NA        NA     NA     NA
#> 8          NA           NA        NA     NA     NA
#> 12         NA           NA        NA     NA     NA
#> 16         NA           NA        NA     NA     NA
#> 20         NA           NA        NA     NA     NA
#> 24         NA           NA        NA     NA     NA
# }
```
