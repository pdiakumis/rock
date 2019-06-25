# Input: Manta VCF
# Output: chrom1, pos1, chrom2, pos2, svtype

#'
#' @examples
#'
#' \dontrun{
#' vcf_file <- "../sv-visualisation/scripts/samplot/samplot/src/2016_249_17_MH_P033-sv-prioritize-manta.vcf.gz"
#' bam_file <- c("bam1", "bam2")
#' sample_nm <- c("s1", "s2")
#' samplot_prep_vcf(vcf_file, bam_file, sample_nm)
#' }
#'
samplot_prep_vcf <- function(vcf_file, bam_file, sample_nm, other = "") {

  stopifnot(length(vcf_file) == 1 && file.exists(vcf_file))
  stopifnot(length(bam_file) == length(sample_nm))
  stopifnot(!any(duplicated(bam_file), duplicated(sample_nm)))
  stopifnot(length(other) == 1)

  prep_manta_vcf(vcf_file)$sv %>%
    dplyr::mutate(
      svlen = ifelse(.data$svtype != "BND", pos2 - pos1, 0),
      cl = paste0(glue::glue("python samplot.py --sv_type {.data$svtype} {other}"),

                  dplyr::case_when(
                    .data$svtype == "BND" ~ glue::glue(
                      "--chrom {.data$chrom1} --start {.data$pos1} --end {.data$pos1} ",
                      "--chrom {.data$chrom2} --start {.data$pos2} --end {.data$pos2} ",
                      "--output_file {.data$svtype}_{.data$chrom1}:{.data$pos1}_{.data$chrom2}:{.data$pos2}"),
                    .data$svtype != "BND" ~ glue::glue(
                      "--chrom {.data$chrom1} --start {.data$pos1} --end {.data$pos2} ",
                      "--output_file {.data$svtype}_{.data$chrom1}:{.data$pos1}-{.data$pos2}"))
      )
    )
}
