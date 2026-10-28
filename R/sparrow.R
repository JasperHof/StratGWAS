#' SPARROW: region- and gene-based (co-)heritability analysis of rare variants
#'
#' Estimates the heritability explained by genes / regions (per variant category)
#' and tests it using a saddlepoint approximation. Phenotypes are matched to the
#' genotypes by IID, covariates are regressed out of the phenotypes here (and out
#' of the genotypes inside the C++ code), and all bookkeeping that used to be
#' manual (cov_df, binary_raw, prs_mask) is done automatically.
#'
#' @param filename PLINK prefix (.bed/.bim/.fam) holding rare variants and common SNPs.
#' @param pheno Data frame or file: FID, IID, then one column per phenotype.
#' @param out_file Output prefix. Writes <out_file> (regions), <out_file>.genes and
#'   <out_file>.windows.
#' @param covariates Optional data frame or file: FID, IID, then numeric covariates
#'   (an intercept is added automatically; do not include one).
#' @param annotation Optional 0/1 matrix or file (with header), one column per variant
#'   category (e.g. rare / common) and one row per .bim variant. A file may start with a
#'   column of variant IDs, which is matched to the .bim. NULL = one category with all variants.
#' @param annot_names Category names (default: column names of `annotation`).
#' @param genes Optional data frame or file: gene, chr, start, end (extra columns ignored).
#'   Must be on the same genome build as the .bim.
#' @param alpha Heritability models tested; alpha[1] is used for the estimates.
#' @param loco_prs Optional LDAK-KVIK Step 1 `.loco.prs` files: one per phenotype
#'   (character vector), or a matrix with one row per phenotype and one column per source.
#' @param kvik_summaries Optional LDAK-KVIK Step 2 `.summaries` files, one per phenotype.
#' @param binary TRUE if the phenotypes are 0/1 (1/2 coding is converted).
#' @param coher TRUE to estimate co-heritability for all pairs of phenotypes.
#' @param h2_without_prs With `loco_prs`: estimate h2 from y (TRUE) and use y - LOCO
#'   only for p-values. Recommended; see the supplement.
#' @param prs_mask Advanced. Categories x PRS sources; which full PRS is subtracted in
#'   the window pass. Default: 0 if h2_without_prs, else 1.
#' @param do_windows Estimate ~1 Mb windows and report window-calibrated h2 (h2_adj).
#' @param chunk_size,window_bp,max_window_chunks Region and window sizes.
#' @param max_gene_snps,min_gene_snps Genes are split into blocks of at most
#'   `max_gene_snps` variants; genes with fewer than `min_gene_snps` are skipped.
#' @param acatv_mac Variants with MAC below this are collapsed into a burden in ACAT-V.
#' @param spa Compute saddlepoint p-values (FALSE = estimates only, much faster).
#' @param chr Optional single chromosome to analyse.
#' @param n_threads Number of threads (0 = RcppParallel default).
#' @param batch_size Units per parallel batch (affects progress reporting only).
#'
#' @return Invisibly, a list with the output file names, the analysed traits and
#'   number of individuals, and the summary returned by the C++ code.
#'
#' @examples
#' \donttest{
#' # example data shipped with the package
#' dir <- system.file("extdata", package = "StratGWAS")
#'
#' res <- sparrow(filename       = file.path(dir, "data"),
#'                pheno          = file.path(dir, "pheno_quant"),
#'                out_file       = file.path(tempdir(), "sparrow_example"),
#'                annotation     = file.path(dir, "annot"),
#'                genes          = file.path(dir, "genes"),
#'                loco_prs       = file.path(dir, "kvik.step1.loco.prs"),
#'                kvik_summaries = file.path(dir, "kvik.step2.summaries"),
#'                n_threads      = 2)
#'
#' genes <- read.table(res$files["genes"], header = TRUE)
#' head(genes[genes$category == "ACAT", ])          # one p-value per gene
#' }
#' @export
sparrow <- function(filename, pheno, out_file,
                    covariates = NULL,
                    annotation = NULL, annot_names = NULL,
                    genes = NULL,
                    alpha = c(-1, -0.5, 0),
                    loco_prs = NULL, kvik_summaries = NULL,
                    binary = FALSE, coher = FALSE,
                    h2_without_prs = TRUE, prs_mask = NULL,
                    do_windows = FALSE, chunk_size = 256, window_bp = 1e6,
                    max_window_chunks = 8,
                    max_gene_snps = 512, min_gene_snps = 2,
                    acatv_mac = 10, spa = TRUE, chr = NULL,
                    n_threads = 0, batch_size = 64) {

  # ---- basic checks -----------------------------------------------------------
  for (ext in c(".bed", ".bim", ".fam"))
    if (!file.exists(paste0(filename, ext))) stop("File not found: ", filename, ext)
  if (!dir.exists(dirname(out_file))) stop("Output directory does not exist: ", dirname(out_file))
  alpha <- unique(as.numeric(alpha))
  if (length(alpha) == 0 || any(!is.finite(alpha))) stop("alpha must contain finite values")

  # ---- phenotypes, matched to the .fam by IID ---------------------------------
  fam <- read.table(paste0(filename, ".fam"), colClasses = "character")
  ph <- .sparrow_read_ids(pheno, "pheno")
  traits <- names(ph)[-(1:2)]
  P <- length(traits)
  if (coher && P < 2) stop("coher = TRUE needs at least two phenotypes")

  keep <- ph$IID %in% fam[[2]]
  message(sum(keep), " of ", nrow(ph), " individuals in `pheno` are in the .fam file")

  # ---- covariates ----------------------------------------------------------------
  C <- NULL
  if (!is.null(covariates)) {
    cv <- .sparrow_read_ids(covariates, "covariates")
    idx <- match(ph$IID, cv$IID)
    C <- as.matrix(cv[idx, -(1:2), drop = FALSE])
    ok <- !is.na(idx) & stats::complete.cases(C)
    message("Dropping ", sum(keep & !ok), " individuals with missing covariates")
    keep <- keep & ok
  }

  Y <- as.matrix(ph[, traits, drop = FALSE])
  # individuals without any phenotype are of no use; with a single phenotype, or
  # with co-heritability, require complete phenotypes (padding biases covariances)
  if (P == 1 || coher) ok <- stats::complete.cases(Y) else ok <- rowSums(!is.na(Y)) > 0
  message("Dropping ", sum(keep & !ok), " individuals with missing phenotypes")
  keep <- keep & ok

  Y <- Y[keep, , drop = FALSE]
  ids <- ph$IID[keep]
  if (!is.null(C)) C <- C[keep, , drop = FALSE]
  n <- nrow(Y)
  if (n < 100) stop("Only ", n, " individuals left after matching -- check the IDs")

  # ---- binary phenotypes: 0/1 coding, keep the raw values -----------------------
  B <- NULL
  if (binary) {
    vals <- unique(as.vector(Y[!is.na(Y)]))
    if (all(vals %in% c(1, 2))) {
      message("Binary phenotypes coded 1/2: converting to 0/1 (2 = case)")
      Y <- Y - 1
    } else if (!all(vals %in% c(0, 1))) stop("binary = TRUE but the phenotypes are not 0/1 (or 1/2)")
    for (t in seq_len(P)) {
      prev <- mean(Y[, t], na.rm = TRUE)
      message(sprintf("  %s: %d cases, prevalence %.4f", traits[t], sum(Y[, t] == 1, na.rm = TRUE), prev))
      if (prev == 0 || prev == 1) stop("Phenotype ", traits[t], " has no cases or no controls")
    }
    B <- Y
  }

  # ---- regress covariates (plus intercept) out of every phenotype ---------------
  Z <- if (is.null(C)) matrix(1, n, 1) else cbind(1, C)
  cov_df <- qr(Z)$rank                           # intercept + linearly independent covariates
  if (cov_df < ncol(Z)) message("Covariates are collinear: using ", cov_df, " degrees of freedom")
  n_pad <- 0
  for (t in seq_len(P)) {
    obs <- !is.na(Y[, t])
    r <- stats::lm.fit(Z[obs, , drop = FALSE], Y[obs, t])$residuals
    Y[, t] <- 0                                  # missing values -> residual mean (padding)
    Y[obs, t] <- r
    n_pad <- n_pad + sum(!obs)
  }
  if (n_pad > 0) message("Padded ", n_pad, " missing phenotype values with the mean (after covariate adjustment)")

  # The phenotypes are NOT standardised: vg is on the phenotype scale and h2 = vg / Var(y),
  # so LDAK-KVIK polygenic scores (on the phenotype scale) can be subtracted directly.
  # Exception: binary co-heritability assumes y = residual / sd (the C++ code recovers the
  # gap between cases and controls from that), so those phenotypes are scaled.
  if (binary && coher) Y <- scale(Y, center = FALSE, scale = apply(Y, 2, stats::sd))

  for (t in seq_len(P))
    if (!(stats::var(Y[, t]) > 0)) stop("Phenotype ", traits[t], " has no variance after covariate adjustment")
  rownames(Y) <- ids
  colnames(Y) <- traits

  binary_raw <- NULL
  if (binary && coher) { binary_raw <- B; rownames(binary_raw) <- ids }

  # ---- variant categories -------------------------------------------------------
  # A file (or data frame) with a header and one 0/1 column per category. If the first
  # column holds variant IDs, rows are matched to the .bim; otherwise they must be in
  # .bim order.
  if (!is.null(annotation)) {
    if (is.character(annotation) && length(annotation) == 1) {
      if (!file.exists(annotation)) stop("annotation file not found: ", annotation)
      annotation <- read.table(annotation, header = TRUE, check.names = FALSE)
    }
    if (is.data.frame(annotation) && !is.numeric(annotation[[1]])) {
      bim_id <- read.table(paste0(filename, ".bim"), colClasses = "character")[[2]]
      idx <- match(bim_id, as.character(annotation[[1]]))
      if (anyNA(idx)) message(sum(is.na(idx)), " variants not in the annotation: assigned to no category")
      annotation <- as.matrix(annotation[idx, -1, drop = FALSE])
      annotation[is.na(annotation)] <- 0
    }
    annotation <- as.matrix(annotation)
    if (anyNA(annotation) || !all(annotation %in% c(0, 1))) stop("annotation must contain only 0 and 1")
    if (any(rowSums(annotation) > 1))
      message("Note: some variants are in several categories; each is assigned to the first")
    storage.mode(annotation) <- "integer"
    if (is.null(annot_names)) annot_names <- colnames(annotation)
    if (is.null(annot_names)) annot_names <- paste0("cat", seq_len(ncol(annotation)))
    if (length(annot_names) != ncol(annotation)) stop("annot_names must have one name per annotation column")
  }
  n_cat <- if (is.null(annotation)) 1 else ncol(annotation)

  # ---- genes ------------------------------------------------------------------------
  if (!is.null(genes)) {
    if (is.character(genes) && length(genes) == 1) {
      if (!file.exists(genes)) stop("genes file not found: ", genes)
      genes <- read.table(genes, header = .sparrow_has_header(genes), colClasses = "character")
    }
    genes <- as.data.frame(genes)
    if (ncol(genes) < 4) stop("genes needs 4 columns: gene, chr, start, end")
    num <- function(x) format(as.numeric(as.character(x)), scientific = FALSE, trim = TRUE)
    genes <- cbind(trimws(as.character(genes[[1]])), trimws(as.character(genes[[2]])),
                   num(genes[[3]]), num(genes[[4]]))
  }

  # ---- LDAK-KVIK inputs -------------------------------------------------------------
  if (!is.null(loco_prs)) {
    loco_prs <- as.matrix(loco_prs)
    if (ncol(loco_prs) == P && nrow(loco_prs) == 1 && P > 1) loco_prs <- t(loco_prs)   # vector given
    if (nrow(loco_prs) != P) stop("loco_prs needs one file per phenotype")
    miss <- loco_prs[!file.exists(loco_prs)]
    if (length(miss)) stop("loco_prs file(s) not found: ", paste(miss, collapse = ", "))
    if (is.null(prs_mask))
      prs_mask <- matrix(if (h2_without_prs) 0L else 1L, n_cat, ncol(loco_prs))
    prs_mask <- as.matrix(prs_mask); storage.mode(prs_mask) <- "integer"
    if (nrow(prs_mask) != n_cat || ncol(prs_mask) != ncol(loco_prs))
      stop("prs_mask must be n_categories x n_prs_sources (", n_cat, " x ", ncol(loco_prs), ")")
  } else if (!is.null(prs_mask)) {
    warning("prs_mask is ignored without loco_prs"); prs_mask <- NULL
  }
  if (!is.null(kvik_summaries)) {
    kvik_summaries <- as.character(kvik_summaries)
    if (length(kvik_summaries) != P) stop("kvik_summaries needs one file per phenotype")
    miss <- kvik_summaries[!file.exists(kvik_summaries)]
    if (length(miss)) stop("kvik_summaries file(s) not found: ", paste(miss, collapse = ", "))
  }

  message(sprintf("SPARROW: %d individuals, %d phenotype(s), %d covariate df, alpha = %s",
                  n, P, cov_df, paste(alpha, collapse = ", ")))

  # ---- run ----------------------------------------------------------------------------
  res <- stratgwas_run(
    filename = filename, pheno_mat = Y,
    annotation = annotation, annot_names = annot_names,
    loco_prs = loco_prs, prs_mask = prs_mask,
    chunk_size = chunk_size, window_bp = window_bp, do_windows = do_windows,
    max_window_chunks = max_window_chunks,
    genes = genes, max_gene_snps = max_gene_snps, min_gene_snps = min_gene_snps,
    alpha = alpha, kvik_summaries = kvik_summaries, acatv_mac = acatv_mac,
    h2_without_prs = h2_without_prs,
    covariates = C, cov_df = cov_df,
    SPA = spa, binary = binary, binary_raw = binary_raw, coher = coher,
    chr = chr, out_file = out_file, batch_size = batch_size, n_threads = n_threads)

  files <- c(regions = out_file,
             genes = if (!is.null(genes)) paste0(out_file, ".genes") else NA,
             windows = if (do_windows) paste0(out_file, ".windows") else NA)
  invisible(list(files = files[!is.na(files)], traits = traits, n = n, info = res))
}


# ---- internal helpers -------------------------------------------------------------------

# Does a whitespace-delimited file have a header? (third field of line 1 not a number)
.sparrow_has_header <- function(path) {
  f <- strsplit(trimws(readLines(path, n = 1)), "[[:space:]]+")[[1]]
  length(f) >= 3 && f[3] != "NA" && is.na(suppressWarnings(as.numeric(f[3])))
}

# Read a FID / IID / values table (data frame or file) and check it.
.sparrow_read_ids <- function(x, what) {
  if (is.character(x) && length(x) == 1) {
    if (!file.exists(x)) stop(what, " file not found: ", x)
    hdr <- .sparrow_has_header(x)
    x <- read.table(x, header = hdr, colClasses = "character",
                    na.strings = c("NA", "-9", ""), check.names = FALSE)
    if (!hdr && ncol(x) >= 3)                    # no header: Pheno1, Pheno2, ... / Cov1, ...
      names(x)[-(1:2)] <- paste0(if (what == "pheno") "Pheno" else "Cov", seq_len(ncol(x) - 2))
  }
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  if (ncol(x) < 3) stop(what, " needs columns FID, IID and at least one value column")
  names(x)[1:2] <- c("FID", "IID")
  x$FID <- trimws(as.character(x$FID)); x$IID <- trimws(as.character(x$IID))
  if (anyDuplicated(x$IID)) stop(what, ": duplicated IIDs (", sum(duplicated(x$IID)), ")")
  for (j in 3:ncol(x)) {
    v <- suppressWarnings(as.numeric(as.character(x[[j]])))
    bad <- is.na(v) & !is.na(x[[j]]) & !(as.character(x[[j]]) %in% c("NA", "-9", ""))
    if (any(bad)) stop(what, " column '", names(x)[j], "' is not numeric",
                       if (what == "covariates") " (code categorical covariates as 0/1 dummies)" else "")
    x[[j]] <- v
  }
  x
}
