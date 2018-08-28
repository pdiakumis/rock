#' Export CNV Segments for IGV Viewing
#'
#' Exports an object of class `cnv` to a
#' [bedgraph](http://genome.ucsc.edu/goldenPath/help/bedgraph.html) file for viewing
#' in [IGV](http://software.broadinstitute.org/software/igv/).
#'
#' @param cnv Object of class `cnv`. See examples.
#' @param out_file Path to write output to. **Needs to have a `bedgraph` suffix**.
#' @param track_name Name of track to appear in IGV (default: "cnv_segs").
#' @param col1 Colour name for amplifications (default: darkgreen)
#' @param col2 Colour name for deletions (default: red)
#' @param autoscale With autoscaling on, IGV automatically adjusts the
#'   plot Y scale to the data range currently in view.
#'   As the user pans and moves, this scaling continually adjusts.
#'   Available options: on, off (default: on).
#'
#' @return Invisible list with the modified `cnv` object and the track header, where:
#'   * `chrom`: chromosome
#'   * `start`: start coordinate
#'   * `end`: end coordinate
#'   * `tot_cn - 2`: total copy number estimate minus 2 (so that CN2 will be the baseline)
#'   * `header`: track information
#'
#' @examples
#' \dontrun{
#' cn_fname <- system.file("extdata", "HCC2218_purple.cnv.tsv", package = "pebbles")
#' cnv <- prep_purple_seg(cn_fname)
#' cnv2igv(cnv, out_file = "~/Desktop/tmp/cnv_segs4igv.bedgraph")
#' }
#'
#' @export
cnv2igv <- function(cnv, out_file = NULL, track_name = "cnv_segs", col1 = "darkgreen", col2 = "red", autoscale = "on") {

  stopifnot(class(cnv) == "cnv")
  stopifnot(!is.null(out_file), grepl(".bedgraph$", out_file))
  stopifnot(length(col1) == 1, length(col2) == 1)
  stopifnot(autoscale %in% c("on", "off"))

  cnv <- cnv$cnv %>%
    dplyr::mutate(tot_cn = .data$tot_cn - 2)

  col1 <- grDevices::col2rgb(col1) %>% paste(collapse = ",")
  col2 <- grDevices::col2rgb(col2) %>% paste(collapse = ",")
  track_head <- glue::glue("track color={col1} altColor={col2} name={track_name} autoScale={autoscale}")
  readr::write_lines(x = track_head, path = out_file)
  readr::write_tsv(x = cnv, path = out_file, append = TRUE)

  invisible(list(header = track_head, cnv = cnv))
}
