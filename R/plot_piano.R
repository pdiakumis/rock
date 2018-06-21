#' Generate Piano Plot
#'
#' Generates a piano plot containing copy number variants from multiple samples
#' or callers.
#'
#' @param cnv_list A named list containing one or more `cnv` objects.
#' The element names will be used for plotting the facets so they need to be unique.
#' See the examples for usage.
#' @return A piano plot with copy number variant segments across
#' multiple facets per sample or caller.
#'
#' @examples
#' \dontrun{
#' cnvkit <- system.file("extdata", "HCC2218_cnvkit-call.cns", package = "pebbles")
#' facets <- system.file("extdata", "HCC2218_facets_cncf.tsv", package = "pebbles")
#' titan <- system.file("extdata", "HCC2218_titan.segs.tsv", package = "pebbles")
#' purple <- system.file("extdata", "HCC2218_purple.cnv.tsv", package = "pebbles")
#'
#' cn_facets <- prep_facets_seg(facets)
#' cn_cnvkit <- prep_cnvkit_seg(cnvkit)
#' cn_titan <- prep_titan_seg(titan)
#' cn_purple <- prep_purple_seg(purple)
#'
#'
#' pdf("/path/to/piano1.pdf", width = 7, height = 7)
#' plot_piano(list(cnvkit = cn_cnvkit, facets = cn_facets, purple = cn_purple, titan = cn_titan))
#' dev.off()
#' }
#' @export
plot_piano <- function(cnv_list) {

  multicnv <- .prep_piano(cnv_list)

  # factor for chromosome plot order
  p <- multicnv %>%
    dplyr::mutate(chrom = factor(.data$chrom, levels = c(1:22, "X", "Y")))

  start <- quo(start)
  end <- quo(end)
  y1 <- quo(y1)
  y2 <- quo(y2)
  col <- quo(col)
  var <- quo(var)
  chrom <- quo(chrom)

    ggplot2::ggplot(p) +
    ggplot2::geom_rect(ggplot2::aes(xmin = !!start, xmax = !!end,
                                    ymin = !!y1, ymax = !!y2, fill = !!col)) +
    ggplot2::geom_hline(yintercept = 2, color = "blue", size = 0.25) +
    ggplot2::scale_fill_identity() +
    ggplot2::facet_grid(ggplot2::vars(var), ggplot2::vars(chrom), scales = "free_x") +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = 'white', colour = "grey"),
      strip.text.y = ggplot2::element_text(size = 11),
      strip.text.x = ggplot2::element_text(size = 11),
      panel.spacing.x = grid::unit(0.1, "lines"))

}


# Prepare data for plotting a piano plot
.prep_piano <- function(cnv_list) {

  stopifnot(rlang::is_list(cnv_list))
  stopifnot(rlang::is_dictionaryish(cnv_list))
  stopifnot(all(purrr::map_lgl(cnv_list, function(x) class(x) == "cnv")))

  multicnv <- purrr::map(cnv_list, "cnv") %>%
    dplyr::bind_rows(.id = "var") %>%
    dplyr::mutate(y1 = .data$tot_cn - 0.5,
                  y2 = .data$tot_cn + 0.5,
                  col = dplyr::case_when(
                    .data$tot_cn == 0 ~ "red",
                    .data$tot_cn == 1 ~ "red",
                    .data$tot_cn == 2 ~ "grey95",
                    .data$tot_cn >= 3 ~ "green",
                    TRUE    ~ "purple"))

  multicnv

}
