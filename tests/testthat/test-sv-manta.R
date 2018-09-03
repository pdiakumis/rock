context("Manta VCF")

manta <- system.file("extdata", "HCC2218_manta.vcf", package = "pebbles")
sv_manta <- read_manta_vcf(manta)
bnd_manta <- sv_manta[sv_manta$svtype == "BND", ]
bnd_manta_extra_row <- rbind(bnd_manta, bnd_manta[nrow(bnd_manta),])

sv_manta2 <- prep_manta_vcf(manta)
sv_manta3 <- prep_manta_vcf2(manta)


test_that("BNDs are proper mates", {
  expect_true(nrow(bnd_manta) %% 2 == 0)
  expect_false(nrow(bnd_manta_extra_row) %% 2 == 0)
  expect_true(.manta_proper_pairs(c("a:1", "b:1", "c:1"), c("a:0", "b:0", "c:0")))
})

test_that("object is of sv class", {
  expect_equal(class(sv_manta2), "sv")
})


test_that("chrom starts with hs", {
  expect_true(all(grepl("hs", sv_manta3$chrom1)))
  expect_true(all(grepl("hs", sv_manta3$chrom2)))
})

test_that("col has colors", {
  expect_true(all(grepl("color", sv_manta3$col)))
})
