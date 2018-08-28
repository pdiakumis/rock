context("IGV files")

cn_fname <- system.file("extdata", "HCC2218_purple.cnv.tsv", package = "pebbles")
cnv <- prep_purple_seg(cn_fname)
igv <- cnv2igv(cnv, out_file = tempfile(fileext = "bedgraph"))

test_that("cnv2igv output has correct dim", {
  expect_equal(class(igv), "list")
  expect_true(length(igv) == 2)
  expect_true(length(igv$header) == 1)
})

test_that("track line is correct", {
  expect_true(igv$header == "track color=0,100,0 altColor=255,0,0 name=cnv_segs autoScale=on")
})

test_that("igv cnv output is correct", {
  expect_true(all(colnames(igv$cnv) == c("chrom", "start", "end", "tot_cn")))
})

test_that("cnv2igv requires bedgraph output", {
  expect_error(cnv2igv(cnv, out_file = "foo"))
})
