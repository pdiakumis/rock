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
    tidyr::separate(col = all_cols, into = c("chrom", "start", "end", "id", "mateid", "svtype", "filter"), sep = "\t") %>%
    dplyr::mutate(rowid = 1:nrow(.), chrom1 = chrom, start1 = start, end1 = start, END = end) %>%
    dplyr::select(rowid, chrom1, start1, end1, END, filter, id, mateid, svtype)

  # BNDs
  # just keep first BND mates (i.e. discard duplicated information)
  # see <https://github.com/Illumina/manta/blob/master/docs/developerGuide/ID.md>
  df_bnd <- DF %>%
    dplyr::filter(svtype == "BND") %>%
    dplyr::bind_cols(., .[match(.$id, .$mateid), c("chrom1", "start1")]) %>%
    dplyr::rename(chrom2 = chrom11, start2 = start11) %>%
    dplyr::mutate(end2 = start2,
                  bndid = substring(id, nchar(id))) %>%
    dplyr::filter(bndid == "1")

  # Other
  df_other <- DF %>%
    dplyr::filter(svtype != "BND") %>%
    dplyr::mutate(chrom2 = chrom1, start2 = END, end2 = END)

  # All together now
  svs <- df_other %>%
    dplyr::bind_rows(df_bnd) %>%
    dplyr::select(rowid, chrom1:end1, chrom2:end2, svtype, id, filter) %>%
    dplyr::arrange(rowid)

  return(svs)
}
