context("Manta VCF")

manta <- system.file("extdata", "HCC2218_manta.vcf", package = "pebbles")
sv_manta <- rock:::read_manta_vcf(manta)
bnd_manta <- sv_manta[sv_manta$svtype == "BND", ]
bnd_manta_extra_row <- rbind(bnd_manta, bnd_manta[nrow(bnd_manta),])


test_that("BNDs are proper mates", {
  expect_true(nrow(bnd_manta) %% 2 == 0)
  expect_false(nrow(bnd_manta_extra_row) %% 2 == 0)
  expect_true(rock:::manta_proper_pairs(c("a:1", "b:1", "c:1"), c("a:0", "b:0", "c:0")))
})

sv_manta2 <- prep_manta_vcf(manta)

test_that("object is of sv class", {
  expect_equal(class(sv_manta2), "sv")
})
