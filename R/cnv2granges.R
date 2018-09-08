#' Convert cnv object to GenomicRanges
#'
#' Converts an object of class `cnv` to `GenomicRanges`
#'
#' @param cnv Object of class `cnv`
#' @return A `GenomicRanges` object with the following columns:
#'   * seqnames: chromosome
#'   * ranges: start-end coordinates
#'   * tot_cn: total copy number estimate
#'
#' @examples
#' cnv_fname <- system.file("extdata", "HCC2218_purple.cnv.tsv", package = "pebbles")
#' cnv <- prep_purple_seg(cnv_fname)
#'
#' cnv2gr(cnv)
#'
#' @export
cnv2gr <- function(cnv) {
  stopifnot(class(cnv) == "cnv")
  GenomicRanges::makeGRangesFromDataFrame(cnv$cnv, keep.extra.columns = TRUE, ignore.strand = TRUE)
}
