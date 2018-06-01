#' Read FACETS CNV Segments
#'
#' Reads the output by the FACETS `emcncf` function and exports the `cncf`
#' CNV segment coordinates.
#'
#' @param facets Path to FACETS `emcncf` text file.
#' @return A dataframe (`tibble`) with the following columns:
#'   * chrom: chromosome
#'   * start: start coordinate
#'   * end: end coordinate
#'   * tot_cn: total copy number estimate
#'
#' @examples
#' cn <- system.file("extdata", "HCC2218_facets_cncf.tsv", package = "pebbles")
#' prep_facets_seg(cn)
#'
#' @export
prep_facets_seg <- function(facets) {

  stopifnot(file.exists(facets))

  cnv <- readr::read_tsv(facets, col_types = "cddddddddddddd") %>%
    dplyr::select(.data$chrom, .data$start, .data$end, .data$tcn.em) %>%
    dplyr::rename(tot_cn = .data$tcn.em)

  structure(list(cnv = cnv), class = "cnv")

}

