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
#' \dontrun{
#' prep_cnvkit_seg("/path/to/sample.cnvkit_call.cns")
#' }
#' @export
prep_cnvkit_seg <- function(cnvkit) {
  stopifnot(file.exists(cnvkit), grepl("call.cns$", cnvkit))

  cnv <- readr::read_tsv(cnvkit, col_types = "ciicdddddddd") %>%
    dplyr::select(.data$chromosome, .data$start, .data$end, .data$cn) %>%
    dplyr::rename(tot_cn = .data$cn,
                  chrom = .data$chromosome)

  structure(list(cnv = cnv), class = "cnv")

}

# for debugging
#cnvkit <- "/Users/pdiakumis/Desktop/projects/umccr/tothill_projects/data/a5/structural/E019-cnvkit-call.cns"
