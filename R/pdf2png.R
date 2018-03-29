#' Transform PDF to PNG format
#'
#' Transforms PDF files to PNG files.
#'
#' You need to make sure each file's suffix is
#' correct i.e. 'pdf' for the input PDF file, and 'png' for the
#' output PNG file. If you don't specify the 'png' argument,
#' by default it will output to the same path and with the same name as 'pdf',
#' but with a 'png' suffix.
#'
#'
#' @param pdf Path to input PDF file.
#' @param png Path to output PNG file.
#' @return Path to output PNG file as a side effect.
#' @examples
#' \dontrun{
#' pdf2png("img1.pdf", "img1.png")
#' pdf2png("img1.pdf")
#' }
#' @export
pdf2png <- function(pdf, png = sub("pdf$", "png", pdf)) {
  stopifnot(grepl(".pdf$", pdf))
  stopifnot(file.exists(pdf))
  stopifnot(grepl(".png$", png))

  magick::image_read_pdf(pdf) %>%
    magick::image_write(png, format = "png")
}
