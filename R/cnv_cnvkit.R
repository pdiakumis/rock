#' Read CNVkit CNV Segments
#'
#' Reads the `cnvkit-call.cns` file output by CNVkit and exports the
#' CNV segment coordinates.
#'
#' @param cnvkit Path to CNVkit `cnvkit-call.cns` file.
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


#' Write CNVkit cnr to wig
#'
#' Reads the `cnvkit.cnr` file output by CNVkit and outputs its
#' `depth` column into a wiggle file for easy viewing in IGV.
#'
#' @param cnr Path to CNVkit's `cnvkit.cnr` file, with `chromosome`,
#'   `start`, `end`, `gene`, `log2`, `depth` and `weight` columns.
#' @param out_file Path to write output to. **Needs to have a `wig` suffix**.
#' @param track_name Name of track to appear in IGV (default: "val").
#' @param col Colour name for the bars (default: lightblue).
#'
#' @return The track header invisibly.
#'
#' @examples
#' \dontrun{
#' cnr <- system.file("extdata", "HCC2218-cnvkit_chr22_1K.cnr", package = "pebbles")
#' cnvkit_depth2igv(cnr, out_file = "HCC2218.cov.wig", track_name = "HCC2218_coverage", col = "orange")
#' }
#'
#' @export
cnvkit_depth2igv <- function(cnr, out_file = NULL, track_name = "val", col = "lightblue") {
  stopifnot(file.exists(cnr))
  stopifnot(!is.null(out_file), grepl(".wig$", out_file))
  stopifnot(length(col) == 1)
  stopifnot(is.character(track_name))

  col <- grDevices::col2rgb(col) %>% paste(collapse = ",")
  track_head <- glue::glue("track type=wiggle_0 name={track_name} color={col} span=250")

  x <- readr::read_tsv(cnr, col_types = "ciicddd") %>%
    dplyr::select(.data$chromosome, .data$start, .data$depth) %>%
    dplyr::mutate(depth = round(.data$depth, 2))

  readr::write_lines(x = track_head, path = out_file)

  for (chrom in c(1:22, "X", "Y")) {
    vs <- glue::glue("variableStep chrom={chrom} span=250")
    readr::write_lines(x = vs, path = out_file, append = TRUE)

    x %>%
      dplyr::filter(.data$chromosome == chrom) %>%
      dplyr::select(.data$start, .data$depth) %>%
      readr::write_tsv(path = out_file, append = TRUE)
  }

}
