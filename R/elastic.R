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

#' Perform randomized double-partitioned Haseman-Elston regression
#'
#' @export
he_reg_part <- function (filename, pheno, snp_cat, cat_names, block_size = 1000, outfile) {

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

  hers <- he_multi_part_enrich(filename, multi_he, snp_cat = snp_cat, cat_names = cat_names, block_size = 1000, outfile)

  return(hers)
}