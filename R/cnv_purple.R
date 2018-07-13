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
#' cn <- system.file("extdata", "HCC2218_purple.cnv.tsv", package = "pebbles")
#' prep_purple_seg(cn)
#'
#' @export
prep_purple_seg <- function(purple) {

  cnv <- read_purple_seg(purple) %>%
    dplyr::select(.data$chrom, .data$start, .data$end, .data$tot_cn) %>%
    dplyr::filter(.data$chrom != "MT") %>%
    dplyr::mutate(tot_cn = round(.data$tot_cn, digits = 0))

  structure(list(cnv = cnv), class = "cnv")
}

# read PURPLE cnv file
read_purple_seg <- function(purple) {

  stopifnot(file.exists(purple))

  cnv <- readr::read_tsv(purple, col_types = "ciididdccc") %>%
    dplyr::rename(chrom = .data$`#chromosome`,
                  tot_cn = .data$copyNumber)

  return(cnv)

}

# read PURPLE fitted file
# cn <- system.file("extdata", "HCC2218_purple.fitted.tsv", package = "pebbles")
# read_purple_fitted(cn)
read_purple_fitted <- function(purple) {
  stopifnot(file.exists(purple))

  cnv <- readr::read_tsv(purple, col_types = "ciicddddddddddddddccdddcd") %>%
    dplyr::rename(chrom = .data$`#chromosome`)

  return(cnv)
}
