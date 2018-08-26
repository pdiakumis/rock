context("FACETS tsv")

facets <- system.file("extdata", "HCC2218_facets_cncf.tsv", package = "pebbles")
cn_facets <- prep_facets_seg(facets)
tab_chrom <- table(cn_facets$cnv$chrom)
gr <- cnv2gr(cn_facets)

test_that("object is of cnv class", {
  expect_equal(class(cn_facets), "cnv")
})

test_that("chrom 23 is X", {
  expect_true("X" %in% names(tab_chrom))
  expect_false("23" %in% names(tab_chrom))
})

test_that("cnv can be converted to gr", {
  expect_s4_class(gr, "GRanges")
})
