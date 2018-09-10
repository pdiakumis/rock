context("Get min chrom")


test_that("Min chrom is min", {
  expect_true(min_chrom("1", "2")  == "1")
  expect_true(min_chrom("22", "X")  == "22")
  expect_true(min_chrom("X", "X")  == "X")
  expect_true(min_chrom("X", "3")  == "3")
  expect_true(min_chrom("X", "Y")  == "X")
})

