#' Print current timestamp
#'
#' See <https://stackoverflow.com/a/10740502/2169986>
#'
#' @export
stamp <- function() {
  cat("[", format(Sys.time(), usetz = TRUE), "]", sep = "")
}

# Get minimum chromosome. Alternative to gtools::mixedsort.
min_chrom <- function(chr1, chr2) {
  if (chr1 == chr2) {
    return(chr1)
  } else if (all(c(chr1, chr2) %in% 1:22)) {
    return(min(as.integer(c(chr1, chr2))))
  } else if (all(c(chr1, chr2) %in% c("X", "Y", "MT"))) {
    return(min(c(chr1, chr2)))
  } else if (chr1 %in% 1:22 & chr2 %in% c("X", "Y", "MT")) {
    return(chr1)
  } else if (chr1 %in% c("X", "Y", "MT") & chr2 %in% 1:22) {
    return(chr2)
  } else if (chr2 %in% 1:22 & chr1 %in% c("X", "Y", "MT")) {
    return(chr2)
  } else if (chr2 %in% c("X", "Y", "MT") & chr1 %in% 1:22) {
    return(chr1)
  } else {
    stop(glue::glue("Something went wrong! chr1 is {chr1}; chr2 is {chr2};"))
  }
}

guess_file_type <- function(file) {
  dplyr::case_when(
    grepl("fastq.gz$", file) ~ "FASTQ",
    grepl("fastq$", file) ~ "FASTQ",
    grepl("fq$", file) ~ "FASTQ",
    grepl("fq.gz$", file) ~ "FASTQ",
    grepl("bam$", file) ~ "BAM",
    grepl("sam$", file) ~ "SAM",
    grepl("vcf$", file) ~ "VCF",
    grepl("vcf.gz$", file) ~ "VCF",
    grepl("txt$", file) ~ "TXT",
    grepl("csv$", file) ~ "CSV",
    TRUE ~ "Other")
}

char2num <- function(x, keep = "0-9") {
  return(as.numeric(gsub(paste0("[^", keep, "]+"), "", x)))
}
