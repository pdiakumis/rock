#' Generate Perl Circos Plot Files
#'
#' Generates required files for a Perl circos plot containing
#' Manta structural variants and CNVkit/FACETS/PURPLE/TitanCNA copy number variants.
#'
#' @param outdir Directory to write the files to.
#' @param manta Path to Manta VCF file.
#' @param cnv Path to copy number call file.
#' @return Generates the required files for a Perl circos plot with copy number
#'   variant segments and structural variant links. Returns these files
#'   invisibly.
#'
#' @examples
#' \dontrun{
#' manta <- system.file("extdata", "HCC2218_manta.vcf", package = "pebbles")
#' cnv <- system.file("extdata", "HCC2218_cnvkit-call.cns", package = "pebbles")
#' outdir <- "~/Desktop/tmp/circos"
#'
#' circos_prep(outdir = outdir, manta = manta, cnv = cnv)
#' }
#' @export
circos_prep <- function(outdir = "circos", manta = NULL, cnv = NULL) {

  stopifnot(!is.null(manta), !is.null(cnv))
  stopifnot(all(file.exists(manta, cnv)))
  dir.create(outdir, recursive = TRUE)


  message(glue::glue("Exporting Manta and CNV circos files to '{outdir}'."))
  # prepare Manta/CNVkit circos files
  manta <- prep_manta_vcf2(manta)
  cnv <- prep_cnv_circos(cnv)

  readr::write_tsv(manta, file.path(outdir, "SAMPLE.link.circos"), col_names = FALSE)
  readr::write_tsv(cnv, file.path(outdir, "SAMPLE.cnv.circos"), col_names = FALSE)


  message(glue::glue("Copying circos templates to '{outdir}'."))
  file.copy(system.file("templates/circos", "circos_simple.conf", package = "pebbles"), outdir, overwrite = TRUE)
  file.copy(system.file("templates/circos", "gaps.txt", package = "pebbles"), outdir, overwrite = TRUE)
  file.copy(system.file("templates/circos", "ideogram.conf", package = "pebbles"), outdir, overwrite = TRUE)

  invisible(list(sv = manta, cnv = cnv))


}

#' Generate Perl Circos Plot
#'
#' Generates a Perl circos plot containing
#' Manta structural variants and CNV calls.
#'
#' @param outdir Directory where the files are located, and where the plot will
#'   be written to.
#' @param name Prefix of plot file (suffix: `_circos_cnvkit_manta.png`).
#' @return Generates a Perl circos plot in the outdir.
#'
#' @examples
#' \dontrun{
#' manta <- system.file("extdata", "HCC2218_manta.vcf", package = "pebbles")
#' cnv <- system.file("extdata", "HCC2218_cnvkit-call.cns", package = "pebbles")
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


#' Generate Perl Circos PURPLE BAF Plot
#'
#' Generates a Perl circos plot containing
#' structural variant links, CNV calls, and BAF values from PURPLE.
#'
#' @param purple_dir Directory where the circos files are located ('purple/circos').
#' @param tumor_name Tumor sample name.
#' @param out_dir Directory where the new circos plot will be written,
#'   along with the input/config files.
#' @return Generates a Perl circos plot in the outdir.
#'
#' @examples
#' \dontrun{
#' purple_dir <- '/path/to/purple/circos'
#' out_dir <- '~/Desktop/tmp/circos_new'
#' tumor_name <- 'sampleA_tumor'
#'
#' plot_circos_baf(purple_dir, out_dir, tumor_name)
#' }
#' @export
plot_circos_baf <- function(purple_dir = NULL, out_dir = "circos_baf", tumor_name = NULL) {

  stopifnot(!is.null(purple_dir), !is.null(tumor_name))
  stopifnot(dir.exists(purple_dir))

  req_files <- file.path(purple_dir, paste0(tumor_name, c(".baf.circos", ".cnv.circos", ".map.circos", ".link.circos")))
  stopifnot(all(file.exists(req_files)))
  target_files <- file.path(out_dir, sub(tumor_name, "SAMPLE", basename(req_files)))

  dir.create(out_dir, recursive = TRUE)
  message(glue::glue("Copying circos templates to '{out_dir}'."))
  file.copy(system.file("templates/circos", "circos_baf.conf", package = "pebbles"), out_dir, overwrite = TRUE)
  file.copy(system.file("templates/circos", "gaps.txt", package = "pebbles"), out_dir, overwrite = TRUE)
  file.copy(system.file("templates/circos", "ideogram.conf", package = "pebbles"), out_dir, overwrite = TRUE)

  message(glue::glue("Copying circos input files to '{out_dir}'."))
  file.copy(from = req_files, to = target_files)
  if (Sys.which("circos") != "") {
    cmd <- glue::glue("circos -nosvg -conf {out_dir}/circos_baf.conf -outputdir {out_dir} -outputfile {tumor_name}_circos_baf.png")
    system(cmd)
  } else {
    stop("Can't find 'circos' in your PATH. Exiting.")
  }
}
