#' Prepare Manta VCF for Circos
#'
#' Prepares a Manta VCF for display in a circos plot. Requires bcftools to be
#' installed in your PATH.
#'
#' @param vcf Path to a bgzipped VCF file
#' @return A data.frame (tibble) with the genomic coordinates
#'   of structural variants detected with Manta.
#'
#' @examples
#' \dontrun{
#' prep_manta_vcf("/path/to/sample.vcf.gz")
#' }
#' @export
prep_manta_vcf <- function(vcf) {
# vcf <- "/Users/pdiakumis/Desktop/projects/umccr/tothill_projects/data/a5/vcf/structural/E019-manta.vcf.gz"

  stopifnot(is.character(vcf),
            file.exists(vcf),
            grepl("vcf.gz$", vcf))

  if (system("bcftools -v", ignore.stdout = TRUE) != 0) {
    stop("Oops, bcftools can't be found! Please install/add to PATH.")
  }

  if (system(paste0("bcftools view -h ",  vcf, " | grep 'BPI_START'"), ignore.stdout = TRUE) == 0) {
    message(stamp(), " BPI has been run on this VCF. Using INFO/BPI_START and INFO/BPI_END for coordinates")
    bcftools_query_command <- "bcftools query -f '%CHROM\t%INFO/BPI_START\t%INFO/BPI_END\t%ID\t%INFO/MATEID\t%INFO/SVTYPE\t%FILTER\n'"
  } else {
    message(stamp(), " BPI has not been run on this VCF. Using POS and INFO/END for coordinates")
    bcftools_query_command <- "bcftools query -f '%CHROM\t%POS\t%INFO/END\t%ID\t%INFO/MATEID\t%INFO/SVTYPE\t%FILTER\n'"
  }


  DF <- system(paste(bcftools_query_command, vcf), intern = TRUE) %>%
    tibble::tibble(all_cols = .) %>%
    tidyr::separate(col = .data$all_cols, into = c("chrom", "start", "end", "id", "mateid", "svtype", "filter"), sep = "\t") %>%
    dplyr::mutate(rowid = 1:nrow(.), chrom1 = .data$chrom, start1 = .data$start, end1 = .data$start, END = .data$end) %>%
    dplyr::select(.data$rowid, .data$chrom1, .data$start1, .data$end1, .data$END, .data$filter, .data$id, .data$mateid, .data$svtype)

  # BNDs
  # just keep first BND mates (i.e. discard duplicated information)
  # see <https://github.com/Illumina/manta/blob/master/docs/developerGuide/ID.md>
  df_bnd <- DF %>%
    dplyr::filter(.data$svtype == "BND") %>%
    dplyr::bind_cols(., .[match(.$id, .$mateid), c("chrom1", "start1")]) %>%
    dplyr::rename(chrom2 = .data$chrom11, start2 = .data$start11) %>%
    dplyr::mutate(end2 = .data$start2,
                  bndid = substring(.data$id, nchar(.data$id))) %>%
    dplyr::filter(.data$bndid == "1")

  # Other
  df_other <- DF %>%
    dplyr::filter(.data$svtype != "BND") %>%
    dplyr::mutate(chrom2 = .data$chrom1, start2 = .data$END, end2 = .data$END)

  # All together now
  svs <- df_other %>%
    dplyr::bind_rows(df_bnd) %>%
    dplyr::select(.data$rowid,
                  .data$chrom1, .data$start1, .data$end1,
                  .data$chrom2, .data$start2, .data$end2, .data$svtype, .data$id, .data$filter) %>%
    dplyr::arrange(.data$rowid)

  return(svs)
}
