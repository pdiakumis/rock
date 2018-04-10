#-- Functions for working with TitanCNA --#


.titan_read_segs <- function(fname) {
  readr::read_tsv(fname, col_names = TRUE, col_types = "cciiidcdiciiiid")
}

# See bcbio/structural/titancna.py
.titan_get_sv_type <- function(x) {
  titan_calls <- c(HOMD  = "DEL", DLOH  = "DEL", NLOH  = "LOH", HET   = "CNV",
                   ALOH  = "DUP", GAIN  = "DUP", ASCNA = "DUP", BCNA  = "DUP",
                   UBCNA = "DUP", OUT   = "CNV")
  stopifnot(all(x %in% names(titan_calls)))
  unname(titan_calls[x])
}

# segs_fname <- "sample1_cluster01.segs.txt"
# segs <- .titan_read_segs(segs_fname)
# .titan_get_sv_type(segs$TITAN_call) %>% table()
