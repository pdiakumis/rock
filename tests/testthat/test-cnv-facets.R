context("FACETS tsv")

facets <- system.file("extdata", "HCC2218_facets_cncf.tsv", package = "pebbles")
cn_facets <- prep_facets_seg(facets)

test_that("object is of cnv class", {
  expect_equal(class(cn_facets), "cnv")
})
