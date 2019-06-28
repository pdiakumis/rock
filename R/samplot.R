
#' Prepare Samplot Commands from VCF
#'
#' Prepares Samplot commands from a structural variant VCF and paths
#' to required BAM files.
#'
#' @param vcf_file Path to the Manta VCF file (any SV VCF should work, really). Can be compressed or not.
#' @param bam_file Path(s) to BAM file(s) (atomic vector).
#' @param sample_nm Name(s) of BAM file(s) to display in plot (atomic vector).
#' @param out_dir Output directory.
#' @param other Other arguments to pass to Samplot (single string).
#' @return A dataframe (`tibble`) with the following:
#'   * chrom1: `CHROM`
#'   * pos1: `POS` | `INFO/BPI_START`
#'   * chrom2: `CHROM` (for mate2 if BND)
#'   * pos2: `INFO/END` | `INFO/BPI_END` (for mate1 if BND)
#'   * svtype: `INFO/SVTYPE`
#'   * svlen: `pos2 - pos2`` if not BND, else 0
#'   * cl: Command line for Samplot
#'
#' @examples
#'
#' \dontrun{
#' vcf_file <- "path/to/vcf1.vcf.gz"
#' bam_file <- c("bam1", "bam2")
#' sample_nm <- c("s1", "s2")
#' samplot_commands(vcf_file, bam_file, sample_nm)
#' }
#'
#' @export
samplot_commands <- function(vcf_file, bam_file, sample_nm, out_dir = NULL, other = "") {

  stopifnot(length(vcf_file) == 1 && file.exists(vcf_file))
  stopifnot(length(bam_file) == length(sample_nm))
  stopifnot(!any(duplicated(bam_file), duplicated(sample_nm)))
  stopifnot(length(out_dir) == 1)
  stopifnot(length(other) == 1)

  if (is.null(out_dir)) {
    out_dir <- paste0(Sys.Date(), "_samplot_results")
  }

  prep_manta_vcf(vcf_file)$sv %>%
    dplyr::mutate(
      svlen = ifelse(.data$svtype != "BND", .data$pos2 - .data$pos1, 0),
      zoom = ifelse((.data$svtype == "BND") | (.data$svtype != "BND" & .data$svlen > 10000), "--zoom 20000 ", ""),
      cl = paste0(glue::glue("python samplot.py --sv_type {.data$svtype} ",
                             "--bams {paste(bam_file, collapse = ' ')} ",
                             "--titles {paste(sample_nm, collapse = ' ')} ",
                             "{zoom}",
                             "{other}"),
                  dplyr::case_when(
                    .data$svtype == "BND" ~ glue::glue(
                      "-c {.data$chrom1} -s {.data$pos1} -e {.data$pos1} ",
                      "-c {.data$chrom2} -s {.data$pos2} -e {.data$pos2} ",
                      "-o {out_dir}/{.data$svtype}_{.data$chrom1}:{.data$pos1}_{.data$chrom2}:{.data$pos2}"),
                    .data$svtype != "BND" ~ glue::glue(
                      "-c {.data$chrom1} -s {.data$pos1} -e {.data$pos2} ",
                      "-o {out_dir}/{.data$svtype}_{.data$chrom1}:{.data$pos1}-{.data$pos2}"),
                    TRUE ~ "XXXX")
      )
    )
}
