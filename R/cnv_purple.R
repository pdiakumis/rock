#' Read PURPLE CNV Segments
#'
#' Reads the `sample-purple.cnv` file output by PURPLE and exports the
#' CNV segment coordinates.
#'
#' @param purple Path to PURPLE `sample-purple.cnv` file.
#' @return A dataframe (`tibble`) with the following columns:
#'   * chrom: chromosome
#'   * start: start coordinate
#'   * end: end coordinate
#'   * tot_cn: total copy number estimate
#'
#' @examples
#' purple <- system.file("extdata", "HCC2218_purple.cnv.tsv", package = "pebbles")
#' prep_purple_seg(purple)
#'
#' @export
prep_purple_seg <- function(purple) {

  stopifnot(file.exists(purple))

  cnv <- readr::read_tsv(purple,
                         col_types = readr::cols_only(
                           `#chromosome` = "c", start = "i", end = "i",
                           copyNumber = "d")) %>%
    dplyr::rename(chrom = .data$`#chromosome`,
                  tot_cn = .data$copyNumber) %>%
    dplyr::filter(.data$chrom != "MT")

  structure(list(cnv = cnv), class = "cnv")
}
