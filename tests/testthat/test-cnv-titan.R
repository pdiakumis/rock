context("TitanCNA tsv")

titan <- system.file("extdata", "HCC2218_titan.segs.tsv", package = "pebbles")
cn_titan <- prep_titan_seg(titan)
cn_titan2 <- read_cnv(titan)
gr <- cnv2gr(cn_titan)

test_that("object is of cnv class", {
  expect_equal(class(cn_titan), "cnv")
  expect_equal(class(cn_titan2), "cnv")
})

test_that("cnv can be converted to gr", {
  expect_s4_class(gr, "GRanges")
})
