#' Generate Piano Plot
#'
#' Generates a piano plot containing copy number variants from multiple samples
#' or callers.
#'
#' @param cnv_list A named list containing one or more `cnv` objects.
#' The element names will be used for plotting the facets so they need to be unique.
#' See the examples for usage.
#' @param chromosomes Character vector of chromosomes to plot.
#' @return A piano plot with copy number variant segments across
#' multiple facets per sample or caller:
#'   - X axis: Chromosome Position
#'   - Y axis: Total Copy Number
#'   - Light green: amplifications (up to CN6)
#'   - Dark green: amplifications (greater than CN6, capped off at CN7)
#'   - Red: deletions
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
#' }
#' @export
plot_piano <- function(cnv_list, chromosomes = c(1:22, "X", "Y")) {

  multicnv <- prep_piano(cnv_list, chromosomes)

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

  multicnv %>%
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
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = 'white', colour = "grey85"),
      strip.background = ggplot2::element_rect(colour = "grey95", fill = "grey90"),
      strip.text.y = ggplot2::element_text(size = 11),
      strip.text.x = ggplot2::element_text(size = 8),
      panel.spacing.x = grid::unit(0.01, "lines")
      )

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
prep_piano <- function(cnv_list, chromosomes = c(1:22, "X", "Y")) {

  stopifnot(rlang::is_list(cnv_list))
  stopifnot(rlang::is_dictionaryish(cnv_list))
  stopifnot(all(purrr::map_lgl(cnv_list, function(x) class(x) == "cnv")))
  stopifnot(all(chromosomes %in% c(1:22, "X", "Y")))

  var_names <- names(cnv_list)

  multicnv <- purrr::map(cnv_list, "cnv") %>%
    dplyr::bind_rows(.id = "var") %>%
    dplyr::filter(.data$tot_cn != 2) %>%
    dplyr::mutate(
      y1 = 2,
      y2 = ifelse(.data$tot_cn > 6, 7, .data$tot_cn),
      fill = dplyr::case_when(.data$tot_cn >= 0 & .data$tot_cn < 2 ~ "red",
                              .data$tot_cn > 2 & .data$tot_cn <= 6 ~ "green",
                              .data$tot_cn > 6                     ~ "green3",
                              TRUE                                 ~ "purple"),
      chrom = readr::parse_factor(.data$chrom, levels = c(1:22, "X", "Y")),
      var = readr::parse_factor(.data$var, levels = var_names)) %>%
    dplyr::filter(.data$chrom %in% chromosomes)

  multicnv

}
