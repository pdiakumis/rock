context("CNVkit cns")

cnvkit <- system.file("extdata", "HCC2218_cnvkit-call.cns", package = "pebbles")
cn_cnvkit <- prep_cnvkit_seg(cnvkit)
gr <- cnv2gr(cn_cnvkit)

test_that("object is of cnv class", {
  expect_equal(class(cn_cnvkit), "cnv")
})

test_that("cnv can be converted to gr", {
  expect_s4_class(gr, "GRanges")
})
