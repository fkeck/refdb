
context("Read external")

# test <- read.csv("tests/testthat/data_fragilaria_ncbi.csv")
# test <- read.csv("tests/testthat/data_baetidae_bold.csv")



test_that("Test export to DADA2 formats", {
  test <- read.csv("data_baetidae_bold.csv")
  test_file <- tempfile(fileext = ".fasta")

  test <- refdb_set_fields_BOLD(test)

  refdb_export_dada2(test, test_file, mode = "taxonomy")

  expect_equal(
    readLines(test_file, warn = FALSE),
    readLines("data_baetidae_bold_dada_taxo.fasta", warn = FALSE)
  )

  refdb_export_dada2(test, test_file, mode = "species")

  expect_equal(
    readLines(test_file, warn = FALSE),
    readLines("data_baetidae_bold_dada_species.fasta", warn = FALSE)
  )

})




test_that("Test export to Mothur formats", {
  test <- read.csv("data_baetidae_bold.csv")
  test_file <- tempfile()

  test <- refdb_set_fields_BOLD(test)

  refdb_export_mothur(test, test_file)

  expect_equal(
    readLines(paste0(test_file, ".fasta"), warn = FALSE),
    readLines("data_baetidae_bold_mothur.fasta", warn = FALSE)
  )
  expect_equal(
    readLines(paste0(test_file, ".txt"), warn = FALSE),
    readLines("data_baetidae_bold_mothur.txt", warn = FALSE)
  )

})
