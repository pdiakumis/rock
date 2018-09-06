#' Read CNVkit CNV Segments
#'
#' Reads the `sample-call.cns` file output by CNVkit and exports the
#' CNV segment coordinates.
#'
#' @param cnvkit Path to CNVkit `sample-call.cns` file.
#' @return A dataframe (`tibble`) with the following columns:
#'   * chrom: chromosome
#'   * start: start coordinate
#'   * end: end coordinate
#'   * tot_cn: total copy number estimate
#'
#' @examples
#' cn <- system.file("extdata", "HCC2218_cnvkit-call.cns", package = "pebbles")
#' prep_cnvkit_seg(cn)
#'
#' @export
prep_cnvkit_seg <- function(cnvkit) {

  stopifnot(file.exists(cnvkit))

  cnv <- readr::read_tsv(cnvkit, col_types = "ciicdddddddd") %>%
    dplyr::select(.data$chromosome, .data$start, .data$end, .data$cn) %>%
    dplyr::rename(tot_cn = .data$cn,
                  chrom = .data$chromosome) %>%
    dplyr::filter(.data$chrom != "MT")

  structure(list(cnv = cnv), class = "cnv")

}

#' Prepare CNVkit Segments for Perl Circos
#'
#' Reads the `sample-call.cns` file output by CNVkit and exports the
#' CNV segment coordinates for plotting in Perl Circos.
#'
#' @param cnvkit Path to CNVkit `sample-call.cns` file.
#' @return A dataframe (`tibble`) with the following columns:
#'   * chrom: homo sapiens chromosome
#'   * start: start coordinate
#'   * end: end coordinate
#'   * value: total copy number estimate, minus 2
#'
#' @examples
#' cn <- system.file("extdata", "HCC2218_cnvkit-call.cns", package = "pebbles")
#' prep_cnvkit_circos(cn)
#'
#' @export
prep_cnvkit_circos <- function(cnvkit) {

  cnv <- prep_cnvkit_seg(cnvkit)$cnv
  cnv %>%
    dplyr::mutate(chrom = paste0("hs", .data$chrom),
                  tot_cn = .data$tot_cn - 2) %>%
    dplyr::rename(value = .data$tot_cn)
}
