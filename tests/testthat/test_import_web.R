#    test <- read.csv("tests/testthat/data_toclean_bold.csv")

test_that("Test download", {
  #skip_on_cran()
  qb <- refdb_import_BOLD("Baetis melanonyx", ncbi_taxo = FALSE)
  qb <- refdb_import_BOLD("Baetis melanonyx", ncbi_taxo = TRUE)
  query <- "Fragilaria capucina rbcl"
  qn <- refdb_import_NCBI(query)

})
