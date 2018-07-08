#' Generate Ideogram Plot
#'
#' Generates an ideogram plot for the specified chromosome
#'
#' @param chrom Chromosome character from 1-22, X or Y.
#' @return An ideogram plot for the specified chromosome.
#'
#' @examples
#' \dontrun{
#' plot_ideogram(chrom = "13")
#' }
#' @export
plot_ideogram <- function(chrom) {

  stopifnot(length(chrom) == 1)
  stopifnot(chrom %in% c(1:22, "X", "Y"))

  cytoband_data <- pebbles::cytoband_data

  colours <- c(
    acen = "#D92F27",
    gneg = "#FFFFFF",
    gpos25 = "#C8C8C8",
    gpos50 = "#C8C8C8",
    gpos75 = "#828282",
    gpos100 = "#000000",
    gvar = "#DCDCDC",
    stalk = "#647FA4"
    )

  cytoband_data$col <- colours[cytoband_data$stain]

  start <- quo(start)
  end <- quo(end)
  col <- quo(col)

  p <- cytoband_data %>%
    dplyr::filter(.data$chrom == paste0("chr", !!chrom)) %>%
    ggplot2::ggplot() +
    ggplot2::geom_rect(ggplot2::aes(xmin = !!start, xmax = !!end,
                                    ymin = -2, ymax = -1, fill = !!col)) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::scale_fill_identity()

  base_theme <- ggplot2::theme_bw() +
    ggplot2::theme(
      axis.ticks = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )

  p + base_theme
}
