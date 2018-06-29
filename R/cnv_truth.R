#' Read HCC2218 CNV Truth Set
#'
#' Reads the `HCC2218_truthset_cnv_bcbio.tsv` file in the `pebbles` package and
#' exports the CNV segment coordinates.
#'
#' @param truth Path to truth set.
#' @return A dataframe (`tibble`) with the following columns:
#'   * chrom: chromosome
#'   * start: start coordinate
#'   * end: end coordinate
#'   * tot_cn: total copy number estimate
#'
#' @examples
#' cn <- system.file("extdata", "HCC2218_truthset_cnv_bcbio.tsv", package = "pebbles")
#' prep_truth_seg(cn)
#'
#' @export
prep_truth_seg <- function(truth) {

  cnv <- readr::read_tsv(truth, col_types = "ciid") %>%
    dplyr::filter(.data$chrom != "MT")

  structure(list(cnv = cnv), class = "cnv")
}
