context("PURPLE tsv")

purple <- system.file("extdata", "HCC2218_purple.cnv.tsv", package = "pebbles")
cn_purple <- prep_purple_seg(purple)

test_that("object is of cnv class", {
  expect_equal(class(cn_purple), "cnv")
})
