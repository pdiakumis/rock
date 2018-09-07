#' Export CNV Segments for IGV Viewing
#'
#' Exports an object of class `cnv` to a
#' [bedgraph](http://genome.ucsc.edu/goldenPath/help/bedgraph.html) file for viewing
#' in [IGV](http://software.broadinstitute.org/software/igv/) as bars.
#'
#' @param cnv Either a BED file with `chrom`, `start`, `end` and `tot_cn` columns, or a `cnv` object.
#' @param out_file Path to write output to. **Needs to have a `bedgraph` suffix**.
#' @param track_name Name of track to appear in IGV (default: "cnv_segs").
#' @param col1 Colour name for amplifications (default: darkgreen).
#' @param col2 Colour name for deletions (default: red).
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

  if (class(cnv) != "cnv") {
    if (!file.exists(cnv)) {
      stop("cnv needs to be either a cnv object or an existing BED file!")
    }
  }
  stopifnot(!is.null(out_file), grepl(".bedgraph$", out_file))
  stopifnot(length(col1) == 1, length(col2) == 1)
  stopifnot(autoscale %in% c("on", "off"))

  if (class(cnv) != "cnv") {
    # cnv is a BED
    # check if there's a header
    cnames <- c("chrom", "start", "end", "tot_cn")
    try_head <- readr::read_lines(cnv, n_max = 1)
    if (grepl("start", try_head)) {
      cnv <- readr::read_tsv(cnv, col_types = "ciid") %>%
        purrr::set_names(cnames)
    } else {
      cnv <- readr::read_tsv(cnv, col_types = "ciid", col_names = cnames)
    }
    cnv <- cnv %>% dplyr::filter(.data$chrom != "MT")
    cnv <- structure(list(cnv = cnv), class = "cnv")

  }

  cnv <- cnv$cnv %>%
    dplyr::mutate(tot_cn = .data$tot_cn - 2,
                  tot_cn = round(.data$tot_cn, 2))


  col1 <- grDevices::col2rgb(col1) %>% paste(collapse = ",")
  col2 <- grDevices::col2rgb(col2) %>% paste(collapse = ",")
  track_head <- glue::glue("track color={col1} altColor={col2} name={track_name} autoScale={autoscale}")
  readr::write_lines(x = track_head, path = out_file)
  readr::write_tsv(x = cnv, path = out_file, append = TRUE)

  invisible(list(header = track_head, cnv = cnv))
}

#' Export BED values for IGV Viewing
#'
#' Given a BED-like file with chromosome, start, end and value columns,
#' export a properly formatted file for viewing the values as a scatter plot
#' in [IGV](http://software.broadinstitute.org/software/igv/). Note that
#' `start` needs to be equal to `end`, since we're looking at specific SNPs,
#' not segments.
#'
#' @param bed BED file with `chrom`, `start`, `end` and `value` columns.
#' @param out_file Path to write output to. **Needs to have an `igv` suffix**.
#' @param track_name Name of track to appear in IGV (default: "val").
#' @param col Colour name for the points (default: blue).
#'
#' @return Invisible list with the modified BED file and the track header, where:
#'   * `chrom`: chromosome
#'   * `start`: start coordinate
#'   * `end`: end coordinate
#'   * `name`: just an ID for the row (`chrom:start-end`)
#'   * `val`: value
#'   * `header`: track information
#'
#' @examples
#' \dontrun{
#' bed <- system.file("extdata", "HCC2218_baf.tsv", package = "pebbles")
#' bedval2igv(bed, out_file = "~/Desktop/tmp/baf1.igv", track_name = "hcc2218_baf", col = "red")
#' }
#'
#' @export
#'
bedval2igv <- function(bed, out_file = NULL, track_name = "val", col = "purple") {
  stopifnot(file.exists(bed))
  stopifnot(!is.null(out_file), grepl(".igv$", out_file))
  stopifnot(length(col) == 1)
  stopifnot(is.character(track_name))

  # check if there's a header
  cnames <- c("chrom", "start", "end", "value")
  try_head <- readr::read_lines(bed, n_max = 1)
  if (grepl("start", try_head)) {
    bed <- readr::read_tsv(bed, col_types = "ciid") %>%
      purrr::set_names(cnames)
  } else {
    bed <- readr::read_tsv(bed, col_types = "ciid", col_names = cnames)
  }


  bed <- bed %>%
    dplyr::mutate(value = round(.data$value, 2),
                  name = glue::glue("chr{chrom}:{start}-{end}")) %>%
    dplyr::select(.data$chrom, .data$start, .data$end, .data$name, .data$value)


  col <- grDevices::col2rgb(col) %>% paste(collapse = ",")
  track_head <- glue::glue("#track type=IGV graphType=points name={track_name} color={col}")
  readr::write_lines(x = track_head, path = out_file)
  readr::write_tsv(x = bed, path = out_file, append = TRUE)

  invisible(list(header = track_head, bed = bed))
}
