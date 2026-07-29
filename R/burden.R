#' Compute burdens
#' 
#' @export
burden <- function(bed_prefix, gene_file = NULL, kb_size = NULL, 
                   n_snps_per_window = NULL, target_mac_per_ind = NULL, 
                   out_file = NULL, write_buffer_genes = 500, effects_file = "", 
                   write_burden_scores = TRUE) {

  if (!is.null(gene_file)) {
    compute_gene_burden(bed_prefix       = bed_prefix,
                        gene_file        = gene_file,
                        out_file         = out_file,
                        write_buffer_genes = 500
    )
  }
  if (!is.null(kb_size)) {
    compute_burden_windows(filename, out_file,
                           kb_size = kb_size,
                           effects_file = effects_file,
                           write_burden_scores = write_burden_scores)
  }
  if (!is.null(n_snps_per_window)) {
    compute_burden_windows(filename, out_file,
                           n_snps_per_window = n_snps_per_window,
                           effects_file = effects_file,
                           write_burden_scores = write_burden_scores)
  }
  if (!is.null(target_mac_per_ind)) {
    compute_burden_windows(filename, out_file,
                           target_mac_per_ind = target_mac_per_ind,
                           effects_file = effects_file,
                           write_burden_scores = write_burden_scores)
  }
}

#' Compute burden correlations
#'
#' @export
burden_correlations <- function(burden_file,
                                max_neighbors = 100,
                                n_inds_hint = -1,
                                return_pairs = FALSE) {

  weights <- compute_burden_weights_blockwise(
    burden_file,
    max_neighbors = max_neighbors
  )

  return(weights)
}

#' To a burden test
#'
#' @export
burden_test <- function(bed_prefix, pheno_mat, he_file,
                        category_file = NULL,
                        trait_name = NULL, out_file = NULL,
                        kb_size = 1, min_enrichment = 0.5,
                        write_buffer_size = 500,
                        chunk_size = 5000) {

  burden_enrich_association(
    bed_prefix          = bed_prefix,
    pheno_mat           = pheno_mat,
    he_file             = he_file,
    category_file       = category_file,
    trait_name          = trait_name,
    out_file            = out_file,
    kb_size             = kb_size,
    min_enrichment      = min_enrichment,
    write_buffer_size   = write_buffer_size,
    chunk_size          = chunk_size
  )

}