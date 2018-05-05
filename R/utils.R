#' Print current timestamp
#'
#' See https://stackoverflow.com/a/10740502/2169986
#'
#' @export
stamp <- function() {
  cat("[", format(Sys.time(), usetz = TRUE), "]", sep = "")
}
