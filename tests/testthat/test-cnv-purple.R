context("PURPLE tsv")

purple <- system.file("extdata", "HCC2218_purple.cnv.tsv", package = "pebbles")
cn_purple <- prep_purple_seg(purple)
cn_purple2 <- read_cnv(purple)

test_that("object is of cnv class", {
  expect_equal(class(cn_purple), "cnv")
  expect_equal(class(cn_purple2), "cnv")
})
