#' Read TitanCNA CNV Segments
#'
#' Reads the `sample-titan.segs.tsv` file output by TitanCNA and exports the
#' CNV segment coordinates.
#'
#' @param titan Path to TitanCNA `sample-titan.segs.tsv` file.
#' @return A dataframe (`tibble`) with the following columns:
#'   * chrom: chromosome
#'   * start: start coordinate
#'   * end: end coordinate
#'   * tot_cn: total copy number estimate
#'
#' @examples
#' cn <- system.file("extdata", "HCC2218_titan.segs.tsv", package = "pebbles")
#' prep_titan_seg(cn)
#'
#' @export
prep_titan_seg <- function(titan) {

  # translate Titan call to SV type
  .titan_get_sv_type <- function(x) {
    titan_calls <- c(HOMD  = "DEL", DLOH  = "DEL", NLOH  = "LOH", HET   = "CNV",
                     ALOH  = "DUP", GAIN  = "DUP", ASCNA = "DUP", BCNA  = "DUP",
                     UBCNA = "DUP", OUT   = "CNV")
    stopifnot(all(x %in% names(titan_calls)))
    unname(titan_calls[x])
  }


  read_titan_seg <- function(titan) {

    stopifnot(file.exists(titan))

    cnv <- readr::read_tsv(titan, col_types = "cciiidcdiciiiid") %>%
      dplyr::select(.data$Chromosome, .data$Start_Position.bp., .data$End_Position.bp., .data$Copy_Number,
                    .data$Median_logR, .data$TITAN_call) %>%
      dplyr::mutate(sv_type = .titan_get_sv_type(.data$TITAN_call)) %>%
      dplyr::rename(chrom = .data$Chromosome,
                    start = .data$Start_Position.bp.,
                    end = .data$End_Position.bp.,
                    tot_cn = .data$Copy_Number)

    return(cnv)
  }


  cnv <- read_titan_seg(titan) %>%
    dplyr::select(.data$chrom, .data$start, .data$end, .data$tot_cn) %>%
    dplyr::filter(.data$chrom != "MT")


  structure(list(cnv = cnv), class = "cnv")

}

