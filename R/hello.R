#' Say Hello to My Little Friend(s)
#'
#' Says hello to all elements in `x`
#'
#' @param x An atomic vector
#' @return A character vector of length equal
#'   to `x`, saying 'Hello' to each of its elements.
#' @examples
#' hello('UMCCR')
#' hello(c('Mollie', 'Walter', 'Dexter'))
#' @export
hello <- function(x) {
  stopifnot(is.atomic(x))
  return(glue::glue("Hello, {x}!"))
}
