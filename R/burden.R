#' Compute burdens
#' 
#' @export
burden <- function(bed_prefix, gene_file = NULL, kb_size = NULL, n_snps_per_window = NULL, target_mac_per_ind = NULL, out_file = NULL, write_buffer_genes = 500) {

  if (!is.null(gene_file)) {
    compute_gene_burden(
        bed_prefix       = bed_prefix,
        gene_file        = gene_file,
        out_file         = out_file,
        write_buffer_genes = 500
    )
  }
  if (!is.null(kb_size)) {
    compute_burden_windows(filename, paste0(outdir, "/burden_10kb"),
                       kb_size = kb_size)
  }
  if (!is.null(n_snps_per_window)) {
    compute_burden_windows(filename, paste0(outdir, "/burden_10kb"),
                       n_snps_per_window = n_snps_per_window)
  }
  if (!is.null(target_mac_per_ind)) {
    compute_burden_windows(filename, paste0(outdir, "/burden_10kb"),
                       target_mac_per_ind = target_mac_per_ind)
  }
}