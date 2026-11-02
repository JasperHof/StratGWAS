#' Fit elastic net model
#'
#' Development version: fits a variational Bayes elastic-net polygenic score.
#'
#' @param filename PLINK prefix (.bed/.bim/.fam) of the genotype data
#' @param pheno_mat Numeric matrix of phenotypes (one column per phenotype), with IIDs as rownames
#' @param h2 Assumed SNP heritability (default: 0.3)
#' @param alpha_param Selection parameter for the dependency of SNP heritability on MAF (default: -0.25)
#' @param p_en,F_en Elastic-net prior hyperparameters: proportion of SNPs in the large-effect
#'   component, and the fraction of heritability they explain (defaults: 0.1 and 0.5)
#' @param chunk_size Number of SNPs read per block (default: 256)
#' @param tol Convergence tolerance (default: 1e-6)
#' @param max_scans Maximum number of passes through the data (default: 10)
#' @return List returned by the C++ elastic-net solver
#' @export
elastic <- function(filename,
                    pheno_mat,
                    h2          = 0.3,
                    alpha_param = -0.25,
                    p_en        = 0.1,
                    F_en        = 0.5,
                    chunk_size  = 256,
                    tol         = 1e-6,
                    max_scans   = 10) {


  # Read in chromosome numbers
  bim <- read.table(paste0(filename, ".bim"))
  chr_all <- bim[, 1]

  # Read in MAFs
  cat("Reading in MAFs... \n")
  maf_all <- compute_maf_all(filename)
  
  cat("Fitting elastic net... \n")
  result <- vb_elastic_net_prs_multi(
    filename    = filename,
    pheno_mat   = pheno_mat,
    maf_all     = maf_all,
    chr_all     = chr_all,
    h2          = h2,
    alpha_param = alpha_param,
    p_en        = p_en,
    F_en        = F_en,
    chunk_size  = chunk_size,
    tol         = tol,
    max_scans   = max_scans,
    loco        = FALSE,
    loco_chr    = 0L
  )

  return(result)
}

#' Compute polygenic scores
#'
#' @param filename PLINK prefix (.bed/.bim/.fam) of the genotype data
#' @param snp_input Character vector of SNP IDs
#' @param beta_input Numeric vector of effect sizes, in the same order as `snp_input`
#' @param block_size Number of SNPs read per block (default: 256)
#' @return Data frame with a polygenic score per individual
#' @export
scores <- function(filename,
                   snp_input,
                   beta_input,
                   block_size = 256) {
  prs <- compute_prs(
    filename   = filename,
    snp_input  = snp_input,
    beta_input = beta_input,
    block_size = block_size
  )

  return(prs)
}

#' SPA solver for HE-based variance components with correction
#'
#' Development version of \code{\link{sparrow}}; prefer \code{sparrow()} for new analyses.
#'
#' @param filename PLINK prefix (.bed/.bim/.fam) holding the variants
#' @param pheno Data frame: FID, IID, then one column per phenotype
#' @param out_file Output prefix
#' @param annotation Optional 0/1 integer matrix, one row per .bim variant and one column per category
#' @param annot_names Optional category names
#' @param loco_prs Optional matrix of LDAK-KVIK Step 1 `.loco.prs` files
#' @param prs_mask Optional categories x PRS sources matrix: which full PRS is subtracted
#' @param chunk_size,window_bp,max_window_chunks Region and window sizes
#' @param do_windows Also estimate ~1 Mb windows (default: TRUE)
#' @param alpha Heritability model(s) tested (default: -1)
#' @param covariates Optional numeric covariate matrix
#' @param cov_df Degrees of freedom used by covariates (NA = determined automatically)
#' @param SPA Compute saddlepoint p-values (default: TRUE)
#' @param spa_pval_threshold Apply SPA only to tests with p-values below this threshold
#' @param binary TRUE if the phenotypes are binary
#' @param binary_raw Optional matrix of the raw 0/1 phenotypes (for binary traits)
#' @param coher TRUE to estimate co-heritability for all pairs of phenotypes
#' @param chr Optional single chromosome to analyse
#' @param batch_size Units per parallel batch
#' @param n_threads Number of threads (0 = RcppParallel default)
#' @param genes Optional gene definitions: gene, chr, start, end
#' @param kvik_summaries Optional LDAK-KVIK Step 2 `.summaries` files
#' @param acatv_mac Variants with MAC below this are collapsed into a burden in ACAT-V
#' @return NULL; results are written to `out_file`
#' @export
he_window_spa <- function (filename, pheno, out_file, annotation = NULL, annot_names = NULL,
                          loco_prs = NULL, prs_mask = NULL, chunk_size = 256,
                          window_bp = 1e6, do_windows = T, alpha = -1,
                          covariates = NULL, cov_df = NA, SPA = T,
                          spa_pval_threshold = 0.1, binary = F,
                          binary_raw = NULL, coher = F, chr = NULL,
                          batch_size = 64, n_threads = 0, max_window_chunks = 8,
                          genes = NULL, kvik_summaries = NULL, acatv_mac = 10
                          ) {

  #storage.mode(snp_cat) <- "integer"          # <- force integer
  #snp_cat <- as.matrix(snp_cat)               # ensure it's a matrix, not df
  #cat_names <- as.character(cat_names)

  # Stadardize
  for (j in 3:ncol(pheno)) pheno[, j] <- as.numeric(scale(pheno[, j]))
  for (j in 3:ncol(pheno)) {
    if (any(is.na(pheno[, j]))) {
      pheno[which(is.na(pheno[, j])), j] <- mean(pheno[, j], na.rm = T)
    }
  }
  # Create multivariate phenotype
  multi <- pheno
  colnames(multi) <- c("FID", "IID", paste0("Pheno", seq_len(ncol(multi) - 2)))

  # HE regression on the phenotype
  multi_he <- as.matrix(multi[, -c(1, 2), drop = FALSE])
  rownames(multi_he) <- multi[, 2]

  stratgwas_run(
    filename = filename,
    pheno_mat = multi_he,
    annotation = annotation,
    annot_names = annot_names,
    loco_prs = loco_prs,
    prs_mask = prs_mask,
    chunk_size = chunk_size,
    window_bp = window_bp,
    do_windows = do_windows,
    alpha = alpha,
    covariates = covariates,
    cov_df = cov_df,
    SPA = SPA,
    spa_pval_threshold = spa_pval_threshold,
    binary = binary,
    binary_raw = binary_raw,
    coher = coher,
    chr = chr,
    out_file = out_file,
    batch_size = batch_size,
    n_threads = n_threads,
    max_window_chunks = max_window_chunks,
    genes = genes,
    kvik_summaries = kvik_summaries,
    acatv_mac = acatv_mac
  )

  return(NULL)
}