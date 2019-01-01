#' Read Manta VCF
#'
#' Read main columns of interest from Manta VCF using bcftools
#'
#' This function uses vcfR (https://github.com/knausb/vcfR) or
#' bcftools (https://samtools.github.io/bcftools/bcftools.html) to read in the VCF file.
#' If you've got a large VCF file, or one with humongous gene annotations,
#' vcfR  will probably choke. A much quicker alternative is bcftools, but you need
#' to make sure it has been installed and that it can be found in the R session
#' via the PATH environmental variable. If you're an RStudio user, you can make
#' sure it recognises the user's PATH by opening the RStudio app via the terminal,
#' or perhaps following the suggestions here: https://stackoverflow.com/questions/31121645.
#'
#'
#' @param vcf Path to Manta VCF file. Can be compressed or not.
#' @return A dataframe (`tibble`) with the following fields from the VCF:
#'   * chrom1: `CHROM`
#'   * pos1: `POS` | `INFO/BPI_START`
#'   * pos2: `INFO/END` | `INFO/BPI_END`
#'   * id: `ID`
#'   * mateid: `INFO/MATEID`
#'   * svtype: `INFO/SVTYPE`
#'   * filter: `FILTER`
#'
#' @examples
#' vcf <- system.file("extdata", "HCC2218_manta.vcf", package = "pebbles")
#' rock:::read_manta_vcf(vcf)
#'
read_manta_vcf <- function(vcf) {

  stopifnot(file.exists(vcf))

  # You have two options: use bcftools or vcfR
  # - If bcftools exists, use it
  # - If bcftools doesn't exist, use vcfR
  if (Sys.which("bcftools") == "") {
    # use vcfR
    vcfr <- vcfR::read.vcfR(vcf, verbose = FALSE)
    DF <- tibble::tibble(chrom1 = as.character(vcfR::getCHROM(vcfr)),
                         pos1 = "dummy1",
                         pos2 = "dummy2",
                         id = vcfR::getID(vcfr),
                         mateid = vcfR::extract.info(vcfr, "MATEID"),
                         svtype = vcfR::extract.info(vcfr, "SVTYPE"),
                         filter = vcfR::getFILTER(vcfr))

    if (any(grepl("BPI_START", vcfr@meta))) {
      # use BPI fields
      DF <- dplyr::mutate(DF,
                          pos1 = as.integer(vcfR::extract.info(vcfr, "BPI_START")),
                          pos2 = as.integer(vcfR::extract.info(vcfr, "BPI_END"))
      )
    } else {
      # use typical fields
      DF <- dplyr::mutate(DF,
                          pos1 = as.integer(vcfR::getPOS(vcfr)),
                          pos2 = as.integer(vcfR::extract.info(vcfr, "END")))
    }

  } else {
    # use bcftools
    if (system(paste0("bcftools view -h ",  vcf, " | grep 'BPI_START'"), ignore.stdout = TRUE) == 0) {
      # use BPI fields
      bcftools_query <- "bcftools query -f '%CHROM\t%INFO/BPI_START\t%INFO/BPI_END\t%ID\t%INFO/MATEID\t%INFO/SVTYPE\t%FILTER\n'"
    } else {
      # use typical fields
      bcftools_query <- "bcftools query -f '%CHROM\t%POS\t%INFO/END\t%ID\t%INFO/MATEID\t%INFO/SVTYPE\t%FILTER\n'"
    }

    DF <- system(paste(bcftools_query, vcf), intern = TRUE) %>%
      tibble::tibble(all_cols = .) %>%
      tidyr::separate(col = .data$all_cols,
                      into = c("chrom1", "pos1", "pos2", "id", "mateid", "svtype", "filter"),
                      sep = "\t", convert = TRUE) %>%
      dplyr::mutate(chrom1 = as.character(.data$chrom1))

  }
  return(DF)
}

#' Prepare Manta VCF for Circos
#'
#' Prepares a Manta VCF for display in a Circos plot.
#'
#' This function uses vcfR (https://github.com/knausb/vcfR) or
#' bcftools (https://samtools.github.io/bcftools/bcftools.html) to read in the VCF file.
#' If you've got a large VCF file, or one with humongous gene annotations,
#' vcfR  will probably choke. A much quicker alternative is bcftools, but you need
#' to make sure it has been installed and that it can be found in the R session
#' via the PATH environmental variable. If you're an RStudio user, you can make
#' sure it recognises the user's PATH by opening the RStudio app via the terminal,
#' or perhaps following the suggestions here: https://stackoverflow.com/questions/31121645.
#'
#' @param vcf Path to Manta VCF file. Can be compressed or not.
#' @param filter_pass Keep only variants annotated with a PASS FILTER? (default: FALSE).
#' @return A dataframe (`tibble`) with the following fields from the VCF:
#'   * chrom1: `CHROM`
#'   * pos1: `POS` | `INFO/BPI_START`
#'   * chrom2: `CHROM` (for mate2 if BND)
#'   * pos2: `INFO/END` | `INFO/BPI_END` (for mate1 if BND)
#'   * svtype: `INFO/SVTYPE`. Used for plotting.
#'
#' @examples
#' vcf <- system.file("extdata", "HCC2218_manta.vcf", package = "pebbles")
#' prep_manta_vcf(vcf)
#'
#' @export
prep_manta_vcf <- function(vcf, filter_pass = FALSE) {

  DF <- read_manta_vcf(vcf)

  # BNDs
  # We have POS of mate2 through INFO/END or INFO/BPI_END, just need CHROM.
  # Just keep BND mate with index 1 (i.e. discard duplicated information)
  # see <https://github.com/Illumina/manta/blob/master/docs/developerGuide/ID.md>
  df_bnd <- DF %>%
    dplyr::filter(.data$svtype == "BND") %>%
    dplyr::bind_cols(., .[match(.$id, .$mateid), "chrom1"]) %>%
    dplyr::rename(chrom2 = .data$chrom11) %>%
    dplyr::mutate(bndid = substring(.data$id, nchar(.data$id))) %>%
    dplyr::filter(.data$bndid == "1") %>%
    dplyr::select(.data$chrom1, .data$pos1, .data$chrom2,
                  .data$pos2, .data$id, .data$mateid, .data$svtype, .data$filter)

  stopifnot(.manta_proper_pairs(df_bnd$id, df_bnd$mateid))

  # Non-BNDs
  df_other <- DF %>%
    dplyr::filter(.data$svtype != "BND") %>%
    dplyr::mutate(chrom2 = .data$chrom1) %>%
    dplyr::select(.data$chrom1, .data$pos1, .data$chrom2, .data$pos2,
                  .data$id, .data$mateid, .data$svtype, .data$filter)

  # All together now
  sv <- df_other %>%
    dplyr::bind_rows(df_bnd)

  if (filter_pass) {
    sv <- sv %>%
      dplyr::filter(.data$filter == "PASS")
  }

  sv <- sv %>%
    dplyr::select(.data$chrom1, .data$pos1,
                  .data$chrom2, .data$pos2, .data$svtype)

  structure(list(sv = sv), class = "sv")
}


# Check if Manta BND mates are properly paired
.manta_proper_pairs <- function(id, mid) {
  ext1 <- substring(id, nchar(id))
  ext2 <- substring(mid, nchar(mid))
  pre1 <- substring(id, 1, nchar(id) - 1)
  pre2 <- substring(mid, 1, nchar(mid) - 1)

  # id should end in 1; mateid in 0
  if (all(ext1 == "1") & all(ext2 == "0") & all(pre1 == pre2)) {
    return(TRUE)
  }
  return(FALSE)
}


#' Prepare Manta VCF for Perl Circos
#'
#' Prepares a Manta VCF for display in a Perl circos plot.
#'
#' @inheritParams prep_manta_vcf
#' @param ... Additional arguments for `prep_manta_vcf`.
#' @return A dataframe (`tibble`) with the following fields from the VCF:
#'   * chrom1: homo sapiens `CHROM`
#'   * pos1: `POS` | `INFO/BPI_START`
#'   * pos1b: `POS` | `INFO/BPI_START`
#'   * chrom2: `CHROM` (for mate2 if BND)
#'   * pos2: `INFO/END` | `INFO/BPI_END` (for mate1 if BND)
#'   * pos2b: `INFO/END` | `INFO/BPI_END` (for mate1 if BND)
#'   * col: link colour
#'
#' @export
#' @examples
#' vcf <- system.file("extdata", "HCC2218_manta.vcf", package = "pebbles")
#' prep_manta_vcf2(vcf)
prep_manta_vcf2 <- function(vcf, ...) {
  sv <- prep_manta_vcf(vcf, ...)$sv

  min_chrom_v <- Vectorize(min_chrom)

  chr_cols <- c(hs1 = '(153,102,0)', hs2 = '(102,102,0)', hs3 = '(153,153,30)',
                hs4 = '(204,0,0)', hs5 = '(255,0,0)', hs6 = '(255,0,204)',
                hs7 = '(255,204,204)', hs8 = '(255,153,0)', hs9 = '(255,204,0)',
                hs10 = '(255,255,0)', hs11 = '(204,255,0)', hs12 = '(0,255,0)',
                hs13 = '(53,128,0)', hs14 = '(0,0,204)', hs15 = '(102,153,255)',
                hs16 = '(153,204,255)', hs17 = '(0,255,255)', hs18 = '(204,255,255)',
                hs19 = '(153,0,204)', hs20 = '(204,51,255)', hs21 = '(204,153,255)',
                hs22 = '(102,102,102)', hsX = '(153,153,153)', hsY = '(204,204,204)')

  stopifnot(length(chr_cols) == 24)

  links_coloured <- sv %>%
    dplyr::mutate(min_chrom = paste0("hs", min_chrom_v(.data$chrom1, .data$chrom2))) %>%
    dplyr::mutate(chrom1 = paste0("hs", .data$chrom1),
                  chrom2 = paste0("hs", .data$chrom2)) %>%
    dplyr::mutate(col = dplyr::case_when(
      svtype == "DEL" ~ '(255,0,0)',
      svtype == "DUP" ~ '(0,255,0)',
      svtype == "INS" ~ '(255,0,255)',
      svtype == "INV" ~ '(255,165,0)',
      svtype == "BND" ~ chr_cols[min_chrom],
      TRUE ~ '0,0,0')) %>%
    dplyr::mutate(pos1b = .data$pos1,
                  pos2b = .data$pos2,
                  col = paste0('color=', col)) %>%
    dplyr::select(.data$chrom1, .data$pos1, .data$pos1b, .data$chrom2, .data$pos2, .data$pos2b, .data$col)

  links_coloured
}
