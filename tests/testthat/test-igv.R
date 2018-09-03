context("IGV files")

cn_fname <- system.file("extdata", "HCC2218_purple.cnv.tsv", package = "pebbles")
cnv <- prep_purple_seg(cn_fname)
igv <- cnv2igv(cnv, out_file = tempfile(fileext = "bedgraph"))

bed <- system.file("extdata", "COLO829_chr21_baf.tsv", package = "pebbles")
igv2 <- bedval2igv(bed, out_file = tempfile(fileext = "igv"), track_name = "colo829_baf", col = "red")

test_that("igv output has correct dim", {
  expect_equal(class(igv), "list")
  expect_true(length(igv) == 2)
  expect_true(length(igv$header) == 1)
  expect_true(length(igv2) == 2)
  expect_true(length(igv2$header) == 1)
})

test_that("track line is correct", {
  expect_true(igv$header == "track color=0,100,0 altColor=255,0,0 name=cnv_segs autoScale=on")
  expect_true(igv2$header == "#track type=IGV graphType=points name=colo829_baf color=255,0,0")
})

test_that("igv output is correct", {
  expect_true(all(colnames(igv$cnv) == c("chrom", "start", "end", "tot_cn")))
  expect_true(all(colnames(igv2$bed) == c("chrom", "start", "end", "name", "value")))
})

test_that("cnv2igv requires bedgraph output", {
  expect_error(cnv2igv(cnv, out_file = "foo"))
})

test_that("bedval2igv requires igv output", {
  expect_error(bedval2igv(bed, out_file = "bar"))
})
