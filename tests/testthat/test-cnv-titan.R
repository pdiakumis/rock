context("TitanCNA tsv")

titan <- system.file("extdata", "HCC2218_titan.segs.tsv", package = "pebbles")
cn_titan <- prep_titan_seg(titan)

test_that("object is of cnv class", {
  expect_equal(class(cn_titan), "cnv")
})
