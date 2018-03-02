hello <- function(x) {
  stopifnot(is.atomic(x))
  return(glue::glue("Hello, {x}!"))
}
