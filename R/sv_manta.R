#' Read Manta VCF
#'
#' Read main columns of interest from Manta VCF using bcftools or bedr
#'
#' Uses bcftools (https://samtools.github.io/bcftools/bcftools.html)
#' or bedr (https://cran.r-project.org/web/packages/bedr/index.html) to read
#' in the VCF file.
#'
#' @param vcf Path to Manta VCF file (`.vcf.gz` or `.vcf`).
#' @return A dataframe (`tibble`) with the following fields from the VCF:
#'   * chrom1: `CHROM` (remove `chr` prefix (if any) for hg38 compatibility)
#'   * pos1: `POS` | `INFO/BPI_START`
#'   * pos2: `INFO/END` | `INFO/BPI_END`
#'   * id: `ID`
#'   * mateid: `INFO/MATEID`
#'   * svtype: `INFO/SVTYPE`
#'   * filter: `FILTER`
#'
#' @examples
#' vcf <- system.file("extdata", "HCC2218_manta.vcf", package = "pebbles")
#' vcf2 <- system.file("extdata", "manta_no_bpi.vcf", package = "pebbles")
#' rock:::read_manta_vcf(vcf)
#' rock:::read_manta_vcf(vcf2)
#'
read_manta_vcf <- function(vcf) {

  stopifnot(file.exists(vcf), length(vcf) == 1)

  # You have two options: use bcftools (first choice) or bedr
  if (Sys.which("bcftools") == "") {
    # use bedr
    x <- bedr::read.vcf(vcf, split.info = TRUE, verbose = FALSE)
    DF <- tibble::tibble(chrom1 = as.character(x$vcf$CHROM),
                         pos1 = "dummy1",
                         pos2 = "dummy2",
                         id = x$vcf$ID,
                         mateid = x$vcf$MATEID,
                         svtype = x$vcf$SVTYPE,
                         filter = x$vcf$FILTER)

    if (any(grepl("BPI_START", x$header$INFO[, "ID"]))) {
      # use BPI fields
      DF <- dplyr::mutate(DF,
                          pos1 = as.integer(x$vcf$BPI_START),
                          pos2 = as.integer(x$vcf$BPI_END))

    } else {
      # use typical fields
      DF <- dplyr::mutate(DF,
                          pos1 = as.integer(x$vcf$POS),
                          pos2 = as.integer(x$vcf$END))
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

  DF %>%
    dplyr::mutate(chrom1 = as.character(sub("chr", "", .data$chrom1)),
                  pos1 = char2num(.data$pos1),
                  pos2 = char2num(.data$pos2)) # sometimes get '.' in pos2 (e.g. from purple)
}

#' Prepare Manta VCF for Circos
#'
#' Prepares a Manta VCF for display in a Circos plot.
#'
#' @param vcf Path to Manta VCF file. Can be compressed or not.
#' @param filter_pass Keep only variants annotated with a PASS FILTER?
#'   (default: FALSE).
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
  # If no BPI, Manta doesn't annotate BNDs with END. Let's grab it from the mate's POS.
  # Sometimes with post-processing a mate might get filtered out.
  # In these cases just filter out the orphan mate.
  #
  # Keep BND mate with index 1 (i.e. discard duplicated information)
  # see <https://github.com/Illumina/manta/blob/master/docs/developerGuide/ID.md>
  df_bnd <- DF %>%
    dplyr::filter(.data$svtype == "BND") %>%
    dplyr::bind_cols(., .[match(.$id, .$mateid), c("chrom1", "pos1")]) %>%
    dplyr::rename(chrom2 = .data$chrom11) %>%
    dplyr::mutate(pos2 = ifelse(is.na(.data$pos2), .data$pos11, .data$pos2),
                  bndid = substring(.data$id, nchar(.data$id)))

  orphan_mates <- df_bnd %>%
    dplyr::filter(.data$chrom2 %in% NA) %>%
    dplyr::mutate(orphan = paste0(.data$chrom1, ":", .data$pos1)) %>%
    dplyr::pull(.data$orphan)

  df_bnd <- df_bnd %>%
    dplyr::filter(!is.na(.data$chrom2)) %>%
    dplyr::filter(.data$bndid == "1") %>%
    dplyr::select(.data$chrom1, .data$pos1, .data$chrom2,
                  .data$pos2, .data$id, .data$mateid, .data$svtype, .data$filter)

  if (length(orphan_mates) > 0) {
    warning(glue::glue("The following {length(orphan_mates)} orphan BND mates are removed:\n",
                       paste(orphan_mates, collapse = "\n")))
  }

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

  chr_cols <- c(
    hs1 =  '(128,125,186)', hs2 =  '(145,142,179)', hs3 =  '(161,159,173)',
    hs4 =  '(179,176,166)', hs5 =  '(196,193,160)', hs6 =  '(213,210,153)',
    hs7 =  '(230,228,147)', hs8 =  '(202,218,138)', hs9 =  '(175,209,129)',
    hs10 = '(147,199,120)', hs11 = '(120,190,111)', hs12 = '(92,180,102)',
    hs13 = '(65,171,93)', hs14 = '(65,166,110)', hs15 = '(65,162,128)',
    hs16 = '(65,158,145)', hs17 = '(65,154,163)', hs18 = '(65,150,180)',
    hs19 = '(66,146,198)', hs20 = '(76,142,196)', hs21 = '(86,139,194)',
    hs22 = '(97,135,192)', hsX =  '(107,132,190)', hsY =  '(117,128,188)'
  )

  stopifnot(length(chr_cols) == 24)

  links_coloured <- sv %>%
    dplyr::mutate(alt_chrom1 = !.data$chrom1 %in% c(1:22, "X", "Y"),
                  alt_chrom2 = !.data$chrom2 %in% c(1:22, "X", "Y"))

  alt_sv <- links_coloured %>%
    dplyr::filter(.data$alt_chrom1 | .data$alt_chrom2) %>%
    dplyr::mutate(
      num = dplyr::row_number(),
      x = paste0(.data$num, ": ", .data$svtype, ", ",
                 .data$chrom1, ":", .data$pos1, " <-> ", .data$chrom2, ":", .data$pos2)) %>%
    dplyr::pull(.data$x)

  if (length(alt_sv > 0)) {
    warning(glue::glue(
      "The following {length(alt_sv)} SVs mapping to alt regions are removed:\n",
      paste(alt_sv, collapse = "\n")))
  }


  links_coloured <- links_coloured %>%
    dplyr::filter(!(.data$alt_chrom1 | .data$alt_chrom2)) %>%
    dplyr::mutate(min_chr = paste0("hs", min_chrom_v(.data$chrom1, .data$chrom2))) %>%
    dplyr::mutate(chrom1 = paste0("hs", .data$chrom1),
                  chrom2 = paste0("hs", .data$chrom2)) %>%
    dplyr::mutate(col = dplyr::case_when(
      svtype == "DEL" ~ '(255,0,0)',
      svtype == "DUP" ~ '(0,255,0)',
      svtype == "INS" ~ '(255,0,255)',
      svtype == "INV" ~ '(255,165,0)',
      svtype == "BND" ~ chr_cols[.data$min_chr],
      TRUE ~ '0,0,0')) %>%
    dplyr::mutate(pos1b = .data$pos1,
                  pos2b = .data$pos2,
                  col = paste0('color=', col)) %>%
    dplyr::select(.data$chrom1, .data$pos1, .data$pos1b, .data$chrom2, .data$pos2, .data$pos2b, .data$col)

  links_coloured
}
