context("Test utils")


test_that("Min chrom is min", {
  expect_true(min_chrom("1", "2")  == "1")
  expect_true(min_chrom("22", "X")  == "22")
  expect_true(min_chrom("X", "X")  == "X")
  expect_true(min_chrom("X", "3")  == "3")
  expect_true(min_chrom("X", "Y")  == "X")
})

test_that("File type is correct", {
  expect_equal(guess_file_type("f1.vcf.gz"), "VCF")
  expect_equal(guess_file_type("f1.vcf"), "VCF")
  expect_equal(guess_file_type("f1.fastq.gz"), "FASTQ")
  expect_equal(guess_file_type("f1.fq.gz"), "FASTQ")
  expect_equal(guess_file_type("f1.fastq"), "FASTQ")
  expect_equal(guess_file_type("f1.fq"), "FASTQ")
  expect_equal(guess_file_type("f1.bam"), "BAM")
  expect_equal(guess_file_type("f1.sam"), "SAM")
  expect_equal(guess_file_type("f1.vcf"), "VCF")
  expect_equal(guess_file_type("f1.vcf.gz"), "VCF")
  expect_equal(guess_file_type("f1.txt"), "TXT")
  expect_equal(guess_file_type("f1.csv"), "CSV")
})
