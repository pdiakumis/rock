context("PURPLE tsv")

purple <- system.file("extdata", "HCC2218_purple.cnv.tsv", package = "pebbles")
cn_purple <- prep_purple_seg(purple)
cn_purple2 <- read_cnv(purple)
gr <- cnv2gr(cn_purple)

test_that("object is of cnv class", {
  expect_equal(class(cn_purple), "cnv")
  expect_equal(class(cn_purple2), "cnv")
})

test_that("cnv can be converted to gr", {
  expect_s4_class(gr, "GRanges")
})
