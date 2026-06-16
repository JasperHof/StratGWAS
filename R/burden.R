#' Compute burdens
#' 
#' @export
burden <- function(bed_prefix, gene_file, out_file = NULL, write_buffer_genes = 500) {

    compute_gene_burden(
        bed_prefix       = bed_prefix,
        gene_file        = gene_file,
        out_file         = out_file,
        write_buffer_genes = 500
    )

  return(NULL)
}