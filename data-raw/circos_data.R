require(devtools)
require(OmicCircos)
require(dplyr)

data("UCSC.hg19.chr", package = "OmicCircos")
UCSC.hg19.chr <- UCSC.hg19.chr %>%
  dplyr::mutate_if(is.factor, as.character) %>%
  dplyr::mutate(chrom = sub("chr", "", .data$chrom))


# Prepare angles + colors
db <- OmicCircos::segAnglePo(seg.dat = UCSC.hg19.chr, seg = unique(UCSC.hg19.chr$chrom))
seg_name <- db[, "seg.name"]
seg_num <- length(seg_name)

circos_data <- list(db = db, seg_name = seg_name, seg_num = seg_num)
devtools::use_data(circos_data)
