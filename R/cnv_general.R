#' Read CNVkit, FACETS, PURPLE or TitanCNA CNV Segments
#'
#' Reads the CNV segment TSV file output by CNVkit, FACETS, PURPLE or TitanCNA.
#'
#' @param cnv Path to CNV file.
#' @return A `cnv` list containing a dataframe (`tibble`) with the following columns:
#'   * chrom: chromosome
#'   * start: start coordinate
#'   * end: end coordinate
#'   * tot_cn: total copy number estimate
#'
#' @examples
#' fn <- list.files(system.file("extdata", package = "pebbles"), full.name = TRUE)
#' read_cnv(fn[1])
#'
#' @export
read_cnv <- function(cnv) {
  caller <- get_caller(cnv)
  stopifnot(!is.null(caller))


  fl <- list(facets = prep_facets_seg,
             cnvkit = prep_cnvkit_seg,
             purple = prep_purple_seg,
             titan  = prep_titan_seg)

  fl[[caller]](cnv)

}

# fn <- list.files(system.file("extdata", package = "pebbles"), full.name = TRUE)
# sapply(fn, get_caller)
get_caller <- function(cnv) {
  stopifnot(file.exists(cnv))
  h <- readr::read_lines(cnv, n_max = 1)
  caller <- NULL

  if (grepl("probes", h)) {
    caller <- "cnvkit"
  } else if (grepl("tcn\\.em", h)) {
    caller <- "facets"
  } else if (grepl("segmentStartSupport", h)) {
    caller <- "purple"
  } else if (grepl("TITAN", h)) {
    caller <- "titan"
  } else {
    message("Unknown caller. Make sure you're reading CNV segments ",
            "from\nCNVkit (cnvkit-call.cns),\nFACETS (facets_segs.tsv),\n",
            "PURPLE (purple.cnv), or\nTITAN (titan_segs.tsv).")
  }

  return(caller)
}
