context("CNVkit cns")

cnvkit <- system.file("extdata", "HCC2218_cnvkit-call.cns", package = "pebbles")
cn_cnvkit <- prep_cnvkit_seg(cnvkit)
cn_cnvkit2 <- read_cnv(cnvkit)

test_that("object is of cnv class", {
  expect_equal(class(cn_cnvkit), "cnv")
  expect_equal(class(cn_cnvkit2), "cnv")
})
