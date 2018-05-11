#' Prepare Manta VCF for Circos
#'
#' Prepares a Manta VCF for display in a circos plot
#'
#' @param vcf A vcfR object or a path to a VCF file
#' @return A data.frame (tibble) with the genomic coordinates
#'   of structural variants detected with Manta.
#'
#' @examples
#' \dontrun{
#' prep_manta_vcf(vcfr_object)
#' prep_manta_vcf("/path/to/vcf")
#' }
#' @export
prep_manta_vcf <- function(vcf) {
  # Takes in VCF filename or vcfR object
  stopifnot(class(vcf) == "vcfR" || file.exists(vcf))

  if (file.exists(vcf)) {
      vcf <- vcfR::read.vcf(vcf)
  }

  # Get CHROM, POS, ID, END, MATEID, SVTYPE, BPI_START etc. from VCF
  # Instead of POS, we can use BPI_START
  # Instead of raw END for non-BNDs, we can use BPI_END
  DF <- vcf@fix[, c("CHROM", "POS", "ID", "FILTER")] %>%
    tibble::as_tibble() %>%
    tibble::rowid_to_column(var = "rowid") %>%
    dplyr::bind_cols(vcfR::extract_info_tidy(vcf)) %>%
    dplyr::mutate(chrom1 = CHROM,
                  start1 = BPI_START,
                  end1 = BPI_START,
                  END = BPI_END) %>%
    dplyr::select(rowid, chrom1, start1, end1, END, FILTER,
                  ID, MATEID, SVTYPE, SVLEN)

  # BNDs
  df_bnd <- DF %>%
    dplyr::filter(SVTYPE == "BND") %>%
    dplyr::bind_cols(., .[match(.$ID, .$MATEID), c("chrom1", "start1")]) %>%
    dplyr::rename(chrom2 = chrom11, start2 = start11) %>%
    dplyr::mutate(end2 = start2)

  # Other
  df_other <- DF %>%
    dplyr::filter(SVTYPE != "BND") %>%
    dplyr::mutate(chrom2 = chrom1, start2 = END, end2 = END)

  # All together now
  df_final <- df_other %>%
    dplyr::bind_rows(df_bnd) %>%
    dplyr::select(rowid, chrom1:end1, chrom2:end2, SVTYPE, ID, SVLEN, FILTER) %>%
    dplyr::arrange(rowid)

  return(df_final)
}

