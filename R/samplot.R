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
      cl = paste0(glue::glue("python samplot.py --sv_type {.data$svtype} ",
                             "--bams {paste(bam_file, collapse = ' ')} ",
                             "--titles {paste(sample_nm, collapse = ' ')} ",
                             "{other}"),
                  dplyr::case_when(
                    .data$svtype == "BND" ~ glue::glue(
                      "-c {.data$chrom1} -s {.data$pos1} -e {.data$pos1} ",
                      "-c {.data$chrom2} -s {.data$pos2} -e {.data$pos2} ",
                      "-o {.data$svtype}_{.data$chrom1}:{.data$pos1}_{.data$chrom2}:{.data$pos2}"),
                    .data$svtype != "BND" ~ glue::glue(
                      "-c {.data$chrom1} -s {.data$pos1} -e {.data$pos2} ",
                      "-o {.data$svtype}_{.data$chrom1}:{.data$pos1}-{.data$pos2}"),
                    TRUE ~ "XXXX")
      )
    )
}
