context("piano plot")

require(pebbles)

chr_len <- pebbles::chr_info[1:24, c("NCBI", "length")]

test_that("chromosomes have correct names", {
  expect_true(all(c(1:22, "X", "Y") %in% chr_len$NCBI))
  expect_true(nrow(chr_len) == 24)
})

