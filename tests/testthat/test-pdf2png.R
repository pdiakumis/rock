require(ggplot2)
context("pdf2png")

tmp_dir <- tempdir()
out_pdf <- file.path(tmp_dir, "plot1.pdf")
out_foo <- file.path(tmp_dir, "plot1.foo")
out_png <- file.path(tmp_dir, "plot1.png")

p <- ggplot(mtcars) +
  geom_point(aes(mpg, hp))
ggsave(out_pdf)

test_that("pdf ends with 'pdf'", {
  expect_error(pdf2png(out_foo))
  expect_error(pdf2png(out_foo, out_png))
})

test_that("path to output png is returned", {
  expect_equal(pdf2png(out_pdf), out_png)
  expect_equal(pdf2png(out_pdf, out_png), out_png)
})
