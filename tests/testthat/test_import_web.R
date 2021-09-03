#    test <- read.csv("tests/testthat/data_toclean_bold.csv")

test_that("Test download", {
  #skip_on_cran()
  qb <- refdb_import_BOLD("Baetis melanonyx", ncbi_taxo = FALSE)
  expect_s3_class(qb, "data.frame")

  qb <- refdb_import_BOLD("Baetis melanonyx", ncbi_taxo = TRUE)
  expect_s3_class(qb, "data.frame")

  query <- "Fragilaria capucina rbcl"
  qn <- refdb_import_NCBI(query)
  expect_s3_class(qn, "data.frame")
})
