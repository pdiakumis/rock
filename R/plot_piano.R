#' Generate Piano Plot
#'
#' Generates a piano plot containing copy number variants from multiple samples
#' or callers.
#'
#' @param cnv_list A named list containing one or more `cnv` objects.
#' The element names will be used for plotting the facets so they need to be unique.
#' See the examples for usage.
#' @param chromosomes Character vector of chromosomes to plot.
#' @param hide_x_lab Boolean. Should the X axis chromosome position labels be hidden or not?
#' @param seg.col Character vector of colours to use for deletions, amplifications,
#' extreme amplifications and everything else. Default values: red, green , green3, blue.
#' @return A piano plot with copy number variant segments across
#' multiple facets per sample or caller:
#'   - X axis: Chromosome Position
#'   - Y axis: Total Copy Number
#'   - Colour 1 (red): deletions
#'   - Colour 2 (light green): amplifications (up to CN6)
#'   - Colour 3 (dark green): amplifications (greater than CN6, capped off at CN7)
#'   - Colour 4 (blue): everything else
#'
#' @examples
#' \dontrun{
#' cnvkit <- system.file("extdata", "HCC2218_cnvkit-call.cns", package = "pebbles")
#' facets <- system.file("extdata", "HCC2218_facets_cncf.tsv", package = "pebbles")
#' titan <- system.file("extdata", "HCC2218_titan.segs.tsv", package = "pebbles")
#' purple <- system.file("extdata", "HCC2218_purple.cnv.tsv", package = "pebbles")
#' truth <- system.file("extdata", "HCC2218_truthset_cnv_bcbio.tsv", package = "pebbles")
#'
#' cn_facets <- prep_facets_seg(facets)
#' cn_cnvkit <- prep_cnvkit_seg(cnvkit)
#' cn_titan <- prep_titan_seg(titan)
#' cn_purple <- prep_purple_seg(purple)
#' cn_truth <- prep_truth_seg(truth)
#' cnv_list <- list(truth = cn_truth, cnvkit = cn_cnvkit, facets = cn_facets,
#'                  purple = cn_purple, titan = cn_titan)
#'
#' plot_piano(cnv_list)
#' plot_piano(cnv_list, chromosomes = c(1:10))
#' plot_piano(cnv_list, seg.col = c("orange", "lightblue", "blue", "pink"))
#' }
#' @export
plot_piano <- function(cnv_list,
                       chromosomes = c(1:22, "X", "Y"),
                       hide_x_lab = TRUE,
                       seg.col = c("red", "green", "green3", "blue")) {

  multicnv <- prep_piano(cnv_list, chromosomes, seg.col = seg.col)

  # fake data to make sure chromosome limits are kept
  chr_len <- pebbles::chr_info[1:24, c("NCBI", "length")]
  length_ranges <- c(rbind(rep(1, 24), chr_len$length))
  chr_limits <- tibble::tibble(chrom = rep(chr_len$NCBI, each = 2),
                           pos = length_ranges,
                           tot_cn = 2) %>%
    dplyr::filter(.data$chrom %in% chromosomes) %>%
    dplyr::mutate(chrom = readr::parse_factor(.data$chrom, levels = chromosomes))


  start <- quo(start)
  end <- quo(end)
  y1 <- quo(y1)
  y2 <- quo(y2)
  col <- quo(col)
  fill <- quo(fill)
  var <- quo(var)
  chrom <- quo(chrom)
  mid <- quo(mid)
  pos <- quo(pos)
  tot_cn <- quo(tot_cn)

  p <- multicnv %>%
    ggplot2::ggplot() +
    ggplot2::geom_rect(ggplot2::aes(xmin = !!start, xmax = !!end,
                                    ymin = !!y1, ymax = !!y2, fill = !!fill)) +
    ggplot2::geom_hline(yintercept = 2, color = "blue", size = 0.25) +
    ggplot2::scale_fill_identity() +
    ggplot2::coord_cartesian(ylim = c(0, 7.5)) +
    ggplot2::geom_blank(data = chr_limits, ggplot2::aes(x = !!pos, y = !!tot_cn)) +
    ggplot2::facet_grid(ggplot2::vars(var), ggplot2::vars(chrom),
                        scales = "free_x", space = "free_x", switch = "y", drop = TRUE) +
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 10), expand = c(0, 0)) +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 4), labels = scales::comma, expand = c(0, 0))

  base_theme <- ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 30, vjust = 1, hjust = 1),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_rect(fill = 'white', colour = "grey85"),
    strip.background = ggplot2::element_rect(colour = "grey95", fill = "grey90"),
    strip.text.y = ggplot2::element_text(size = 11),
    strip.text.x = ggplot2::element_text(size = 8),
    panel.spacing.x = grid::unit(0.01, "lines")
  )

  hide_x_lab_theme <- ggplot2::theme(
    axis.ticks.x = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_blank())

  if (nrow(multicnv) == 0) {
    text <- "No CNVs to display here!"
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 4, y = 10, size = 8, label = text) +
      ggplot2::theme_bw() +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank(),
                     axis.text = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(),
                     axis.title = ggplot2::element_blank())

      return(p)
  }

  if (hide_x_lab) {
    return(p + base_theme + hide_x_lab_theme)
  } else {
    return(p + base_theme)
  }

}


# Prepare data for plotting a piano plot
#
# Prepares data for plotting a piano plot
#
# env <- new.env()
# env$cnvkit <- system.file("extdata", "HCC2218_cnvkit-call.cns", package = "pebbles")
# env$facets <- system.file("extdata", "HCC2218_facets_cncf.tsv", package = "pebbles")
# env$titan <- system.file("extdata", "HCC2218_titan.segs.tsv", package = "pebbles")
# env$purple <- system.file("extdata", "HCC2218_purple.cnv.tsv", package = "pebbles")
# env$truth <- system.file("extdata", "HCC2218_truthset_cnv_bcbio.tsv", package = "pebbles")
#
# env$cn_facets <- prep_facets_seg(env$facets)
# env$cn_cnvkit <- prep_cnvkit_seg(env$cnvkit)
# env$cn_titan <- prep_titan_seg(env$titan)
# env$cn_purple <- prep_purple_seg(env$purple)
# env$cn_truth <- prep_truth_seg(env$truth)
# cnv_list <- list(truth = env$cn_truth, cnvkit = env$cn_cnvkit, facets = env$cn_facets, purple = env$cn_purple, titan = env$cn_titan)
prep_piano <- function(cnv_list, chromosomes = c(1:22, "X", "Y"), seg.col = c("red", "green", "green3", "blue")) {

  stopifnot(rlang::is_list(cnv_list))
  stopifnot(rlang::is_dictionaryish(cnv_list))
  stopifnot(all(purrr::map_lgl(cnv_list, function(x) class(x) == "cnv")))
  stopifnot(all(chromosomes %in% c(1:22, "X", "Y")))
  stopifnot(length(seg.col) == 4)

  var_names <- names(cnv_list)
  seg.col <- seg.col %>%
    rlang::set_names(c("del", "amp", "ampx", "def"))

  multicnv <- purrr::map(cnv_list, "cnv") %>%
    dplyr::bind_rows(.id = "var") %>%
    dplyr::filter(.data$tot_cn != 2) %>%
    dplyr::mutate(
      y1 = 2,
      y2 = ifelse(.data$tot_cn > 6, 7, .data$tot_cn),
      fill = dplyr::case_when(.data$tot_cn >= 0 & .data$tot_cn < 2 ~ seg.col["del"],
                              .data$tot_cn > 2 & .data$tot_cn <= 6 ~ seg.col["amp"],
                              .data$tot_cn > 6                     ~ seg.col["ampx"],
                              TRUE                                 ~ seg.col["def"]),
      chrom = readr::parse_factor(.data$chrom, levels = c(1:22, "X", "Y")),
      var = readr::parse_factor(.data$var, levels = var_names)) %>%
    dplyr::filter(.data$chrom %in% chromosomes)

  multicnv

}
