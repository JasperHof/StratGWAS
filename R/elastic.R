#' Fit elastic net model
#'
#' Development version
#'
#' @export
elastic <- function(filename    = filename,
                    pheno_mat   = pheno_mat,
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
#' @export
scores <- function(filename   = filename,
                   snp_input  = snp,
                   beta_input = beta,
                   block_size = 256) {
  prs <- compute_prs(
    filename   = filename,
    snp_input  = snp_input,
    beta_input = beta_input,
    block_size = 256
  )

  return(prs)
}

#' Perform randomized partitioned Haseman-Elston regression
#' 
#' @export
he_reg <- function(pheno, filename) {

  # Stadardize
  for (j in 3:ncol(pheno)) pheno[, j] <- as.numeric(scale(pheno[, j]))
  for (j in 3:ncol(pheno)) {
    if (any(is.na(pheno[, j]))) {
      pheno[which(is.na(pheno[, j])), j] <- mean(pheno[, j], na.rm = T)
    }
  }
  # Create multivariate phenotype
  multi <- pheno
  colnames(multi) <- c("FID", "IID", "Pheno")

  # HE regression on the phenotype
  multi_he <- as.matrix(multi[, -c(1, 2), drop = FALSE])
  rownames(multi_he) <- multi[, 2]

  her <- he_multi_part(filename, multi_he, 1000)

  return(her)
}

#' Perform randomized double-partitioned window-based Haseman-Elston regression for enrichment analysis AND alpha estimation
#' Alpha from -1, -.75, -.5, -.25, 0, by default (and can't change in current implementation)
#'
#' @export
he_reg_part_alpha <- function (filename, pheno, snp_cat, cat_names, common_filename = NULL, alpha = -1,
                               window_size = 1e6, nmcmc = 20, outfile, n_threads = 1,
                               batch_size = 16, genome_wide = FALSE, min_snps = 1000,
                               weights = NULL, covariates = NULL, coher = FALSE,
                               se = FALSE, alpha_common = -1, collapse_mac = 0,
                               collapse_n = 5) {


  storage.mode(snp_cat) <- "integer"          # <- force integer
  snp_cat <- as.matrix(snp_cat)               # ensure it's a matrix, not df
  cat_names <- as.character(cat_names)

  # Stadardize
  for (j in 3:ncol(pheno)) pheno[, j] <- as.numeric(scale(pheno[, j]))
  for (j in 3:ncol(pheno)) {
    if (any(is.na(pheno[, j]))) {
      pheno[which(is.na(pheno[, j])), j] <- mean(pheno[, j], na.rm = T)
    }
  }
  # Create multivariate phenotype
  multi <- pheno
  colnames(multi) <- c("FID", "IID", "Pheno")

  # HE regression on the phenotype
  multi_he <- as.matrix(multi[, -c(1, 2), drop = FALSE])
  rownames(multi_he) <- multi[, 2]

  #hers <- he_multi_part_enrich(filename, multi_he, snp_cat = snp_cat, cat_names = cat_names, block_size = 1000, outfile)
  #hers <- he_multi_part_enrich_alpha(filename, multi_he, snp_cat = snp_cat, cat_names = cat_names, block_size = 1000, outfile)
  
  res <- he_sliding_window_part(
    filename    = filename,               # PLINK prefix (WES)
    pheno_mat   = multi_he,
    snp_cat     = snp_cat,
    cat_names   = cat_names,
    window_size = window_size,
    alpha       = alpha,                  # heritability model; -1 = unweighted (default)
    alpha_common = alpha_common,
    common_filename = common_filename,    # optional common-SNP background GRM
    nmcmc       = nmcmc,
    se          = se,
    out_file    = outfile,
    n_threads   = n_threads,
    batch_size  = batch_size,
    seed        = 12345,
    min_snps    = min_snps,
    weights     = weights,
    covariates  = covariates,
    coher       = coher,
    collapse_mac = collapse_mac,
    collapse_n  = collapse_n
  )
  
  # hers$genome_h2 contains the genome estimates per alpha
  # hers$$best_alpha contains best alpha based on Frobenius norm
  
  #info <- cbind(rownames(hers$genome_h2), hers$genome_h2, hers$genome_se, hers$fit_score, hers$best_alpha)
  #rownames(info) <- colnames(info) <- NULL

  return(res)
}

#' SPA solved for randomized double-partitioned window-based Haseman-Elston regression
#'
#' @export
he_reg_spa <- function (filename, pheno, chunk_size = 256, flank_chunks = 1,
                        min_chunk_snps = 5, alpha = -1, alpha_common = -1,
                        common_filename = NULL,
                        common_window = NA, max_common_snps = 256,
                        outfile, batch_size = 4, n_threads = 1,
                        covariates = NULL, SPA = T, spa_pval_threshold = 0.1,
                        annotation = NULL, annot_names = NULL) {

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
  colnames(multi) <- c("FID", "IID", "Pheno")

  # HE regression on the phenotype
  multi_he <- as.matrix(multi[, -c(1, 2), drop = FALSE])
  rownames(multi_he) <- multi[, 2]

  #hers <- he_multi_part_enrich(filename, multi_he, snp_cat = snp_cat, cat_names = cat_names, block_size = 1000, outfile)
  #hers <- he_multi_part_enrich_alpha(filename, multi_he, snp_cat = snp_cat, cat_names = cat_names, block_size = 1000, outfile)

  he_chunk_spa(
    filename    =  filename,               # PLINK prefix (WES)
    pheno_mat   = multi_he,
    chunk_size  = chunk_size,
    flank_chunks = flank_chunks,
    min_chunk_snps = min_chunk_snps,
    alpha       = alpha,
    alpha_common = alpha_common,
    common_filename = common_filename,
    common_window = common_window,
    max_common_snps = max_common_snps,
    out_file    = outfile,
    batch_size  = batch_size,
    n_threads   = n_threads,
    covariates  = covariates,
    SPA         = T,
    spa_pval_threshold = spa_pval_threshold,
    annotation = annotation,
    annot_names = annot_names
  )

  # hers$genome_h2 contains the genome estimates per alpha
  # hers$$best_alpha contains best alpha based on Frobenius norm
  
  #info <- cbind(rownames(hers$genome_h2), hers$genome_h2, hers$genome_se, hers$fit_score, hers$best_alpha)
  #rownames(info) <- colnames(info) <- NULL

  return(NULL)
}
