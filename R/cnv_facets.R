#' Read FACETS CNV Segments
#'
#' Reads the output by the FACETS `emcncf` function and exports the `cncf`
#' CNV segment coordinates.
#'
#' @param facets Path to FACETS `emcncf` object.
#' @return A dataframe (`tibble`) with the following columns:
#'   * chrom: chromosome
#'   * start: start coordinate
#'   * end: end coordinate
#'   * tot_cn: total copy number estimate
#'
#' @examples
#' \dontrun{
#' prep_facets_seg("/path/to/sample.facets_emcncf.rds")
#' }
#' @export
prep_facets_seg <- function(facets) {
  stopifnot(file.exists(facets), grepl("rds$", facets))

  cnv <- readr::read_rds(facets)[["cncf"]] %>%
    dplyr::select(.data$chrom, .data$start, .data$end, .data$tcn.em) %>%
    dplyr::rename(tot_cn = .data$tcn.em) %>%
    dplyr::mutate(chrom = as.character(.data$chrom)) %>%
    tibble::as_tibble()

  structure(list(cnv = cnv), class = "cnv")

}

# for debugging
# facets <- "/Users/pdiakumis/Desktop/projects/umccr/tothill_projects/data/a5/facets/batch1/E019/E019_cval_150_fit.rds"

