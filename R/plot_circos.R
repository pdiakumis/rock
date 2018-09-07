#' Generate Circos Plot
#'
#' Generates an OmicCircos circos plot containing structural and copy number variants.
#'
#' See the example for how to generate `sv` and `cnv` objects.
#'
#' @param sv An `sv` object
#' @param cnv A `cnv` object
#' @return An OmicCircos circos plot with copy number variant segments
#'   and structural variant links.
#'
#' @examples
#' \dontrun{
#' manta_vcf <- system.file("extdata", "HCC2218_manta.vcf", package = "pebbles")
#' cnvkit_seg <- system.file("extdata", "HCC2218_cnvkit-call.cns", package = "pebbles")
#'
#' sv <- prep_manta_vcf(manta_vcf)
#' cnv <- prep_cnvkit_seg(cnvkit_seg)
#'
#' pdf("/path/to/circos1.pdf", width = 7, height = 7)
#' plot_circos(sv = sv, cnv = cnv)
#' dev.off()
#' }
#' @export
plot_circos <- function(sv = NULL, cnv = NULL) {
  # at least one is required
  stopifnot((!is.null(sv) && class(sv) == "sv") ||
              (!is.null(cnv) && class(cnv) == "cnv"))

  # Prepare angles + colors
  chr_colors <- grDevices::rainbow(pebbles::circos_data$seg_num, alpha = 0.5) %>%
    rlang::set_names(pebbles::circos_data$seg_name)

  #---- SV data ----#
  # var1, var2 just dummy vars for OmicCircos
  if (!is.null(sv)) {

    sv_all <- sv$sv %>%
      dplyr::mutate(var1 = "var1", var2 = "var2") %>%
      dplyr::select(.data$chrom1, .data$pos1, .data$var1,
                    .data$chrom2, .data$pos2, .data$var2,
                    .data$svtype) %>%
      dplyr::mutate(col = dplyr::case_when(
        .data$svtype == "DEL" ~ "red",
        .data$svtype == "DUP" ~ "green",
        .data$svtype == "INS" ~ "purple",
        .data$svtype == "INV" ~ "orange",
        TRUE                  ~ chr_colors[.data$chrom1]))

    # OmicCircos doesn't like tibbles...
    svs_bnd <- sv_all %>% dplyr::filter(.data$svtype == "BND") %>% as.data.frame()
    svs_other <- sv_all %>% dplyr::filter(.data$svtype != "BND") %>% as.data.frame()
  }

  #---- CNV data ----#
  if (!is.null(cnv)) {

    cnv_all <- cnv$cnv %>%
      dplyr::mutate(
        col = dplyr::case_when(
          tot_cn == 0 ~ "red",
          tot_cn == 1 ~ "red",
          tot_cn == 2 ~ "grey50",
          tot_cn >= 3 ~ "green",
          TRUE    ~ "black")) %>%
      as.data.frame()
  } else {
    cnv_all <- data.frame() # need to keep the logic happy
  }

  #---- Circos Plot ----#
  # turn off warnings because in their code they take max/min of character matrix...
  # options(warn = -1)
  graphics::par(mar = c(.5, .5, .5, .5))
  graphics::plot(c(1, 800), c(1, 800), type = "n", axes = FALSE, xlab = "", ylab = "", main = "")

  OmicCircos::circos(R = 400, cir = pebbles::circos_data$db, type = "chr",
                     col = chr_colors, print.chr.lab = TRUE, W = 4)

  max_cnv <- 10000
  if (!is.null(cnv) && nrow(cnv_all) < max_cnv) {
    OmicCircos::circos(R = 260, cir = pebbles::circos_data$db, type = "arc",
                       W = 120, mapping = cnv_all, col.v = 4, B = TRUE, lwd = 5,
                       col = cnv_all$col, scale = FALSE)
  }

  # Case if given both SVs and CNVs, and CNVs are few
  if (!is.null(sv) && !is.null(cnv) && nrow(cnv_all) < max_cnv) {
    # some cases where no BNDs or other SVs PASS.
    if (nrow(svs_bnd) > 0) {
      OmicCircos::circos(R = 260, cir = pebbles::circos_data$db, type = "link",
                         W = 40, mapping = svs_bnd, lwd = 2,
                         col = svs_bnd$col)
    }
    if (nrow(svs_other) > 0) {
      OmicCircos::circos(R = 260, cir = pebbles::circos_data$db, type = "link2",
                         W = 20, mapping = svs_other, lwd = 1,
                         col = svs_other$col)
    }
  }

  # Case if given only SVs, or if given both but CNVs are many
  if ((!is.null(sv) && is.null(cnv)) || (!is.null(sv) && !is.null(cnv) && nrow(cnv_all) >= max_cnv)) {

    if (nrow(cnv_all) >= max_cnv) {
      message(glue::glue("Not plotting copy number segments since ",
                         "{nrow(cnv_all)} > 10,000 and OmicCircos chokes!"))
    }

    # some cases where no BNDs or other SVs PASS.
    if (nrow(svs_bnd) > 0) {
      OmicCircos::circos(R = 380, cir = pebbles::circos_data$db, type = "link",
                         W = 40, mapping = svs_bnd, lwd = 2,
                         col = svs_bnd$col)
    }
    if (nrow(svs_other) > 0) {
      OmicCircos::circos(R = 380, cir = pebbles::circos_data$db, type = "link2",
                         W = 20, mapping = svs_other, lwd = 1,
                         col = svs_other$col)
    }
  }

}

#' Generate Perl Circos Plot Files
#'
#' Generates required files for a Perl circos plot containing
#' Manta structural variants and CNVkit copy number variants.
#'
#' @param outdir Directory to write the files to.
#' @param manta Path to Manta VCF file.
#' @param cnvkit Path to CNVkit `cnvkit-call.cns` file.
#' @return Generates the required files for a Perl circos plot with copy number
#'   variant segments and structural variant links.
#'
#' @examples
#' \dontrun{
#' manta <- system.file("extdata", "HCC2218_manta.vcf", package = "pebbles")
#' cnvkit <- system.file("extdata", "HCC2218_cnvkit-call.cns", package = "pebbles")
#' outdir <- "circos"
#'
#' circos_prep(outdir, manta, cnvkit)
#' }
#' @export
circos_prep <- function(outdir = "circos", manta = NULL, cnvkit = NULL) {

  stopifnot(!is.null(manta), !is.null(cnvkit))
  stopifnot(all(file.exists(manta, cnvkit)))
  dir.create(outdir, recursive = TRUE)


  message(glue::glue("Exporting Manta and CNVkit circos files to '{outdir}'."))
  # prepare Manta/CNVkit circos files
  manta <- prep_manta_vcf2(manta)
  cnvkit <- prep_cnvkit_circos(cnvkit)

  readr::write_tsv(manta, file.path(outdir, "SAMPLE.link.circos"), col_names = FALSE)
  readr::write_tsv(cnvkit, file.path(outdir, "SAMPLE.cnv.circos"), col_names = FALSE)


  message(glue::glue("Copying circos templates to '{outdir}'."))
  file.copy(system.file("templates/circos", "circos_simple.conf", package = "rock"), outdir, overwrite = TRUE)
  file.copy(system.file("templates/circos", "gaps.txt", package = "rock"), outdir, overwrite = TRUE)
  file.copy(system.file("templates/circos", "ideogram.conf", package = "rock"), outdir, overwrite = TRUE)


}

#' Generate Perl Circos Plot
#'
#' Generates a Perl circos plot containing
#' Manta structural variants and CNVkit copy number variants.
#'
#' @param outdir Directory where the files are located, and where the plot will
#'   will be written to.
#' @param name Prefix of plot file (suffix: `_circos_cnvkit_manta.png`).
#' @return Generates a Perl circos plot in the outdir.
#'
#' @examples
#' \dontrun{
#' manta <- system.file("extdata", "HCC2218_manta.vcf", package = "pebbles")
#' cnvkit <- system.file("extdata", "HCC2218_cnvkit-call.cns", package = "pebbles")
#' outdir <- "circos"
#'
#' circos_prep(outdir, manta, cnvkit)
#' plot_circos2(outdir, name = "foo")
#' }
#' @export
plot_circos2 <- function(outdir = "circos", name = "x") {

  if (Sys.which("circos") != "") {
    cmd <- glue::glue("circos -nosvg -conf {outdir}/circos_simple.conf -outputdir {outdir} -outputfile {name}_circos_cnvkit_manta.png")
    system(cmd)
  } else {
    stop("Can't find 'circos' in your PATH. Exiting.")
  }
}
