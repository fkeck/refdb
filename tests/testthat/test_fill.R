
# test <- read.csv("tests/testthat/data_baetidae_bold.csv")



test_that("Test fill upstream", {
  test <- read.csv("data_baetidae_bold.csv")
  test_file <- tempfile(fileext = ".fasta")

  test <- refdb_set_fields_BOLD(test)
  test[2, "subfamily_name"] <- NA
  test[2, "family_name"] <- NA
  test[3, "order_name"] <- NA
  test[3, "species_name"] <- NA

  res <- refdb_fill_tax_upstream(test)

  expect_equal(
    res$family_name[2], paste(test$genus_name[2], "undef. family")
  )
  expect_equal(
    res$subfamily_name[2], paste(test$genus_name[2], "undef. subfamily")
  )
  expect_true(is.na(res$species_name[3]))
  expect_false(is.na(res$order_name[3]))

})


test_that("Test fill downstream", {
  test <- read.csv("data_baetidae_bold.csv")
  test_file <- tempfile(fileext = ".fasta")

  test <- refdb_set_fields_BOLD(test)
  test[2, "subfamily_name"] <- NA
  test[2, "family_name"] <- NA
  test[3, "genus_name"] <- NA
  test[3, "species_name"] <- NA
  test[3, "order_name"] <- NA

  res <- refdb_fill_tax_downstream(test)

  expect_equal(
    res$genus_name[3], paste(test$subfamily_name[3], "indet. genus")
  )
  expect_equal(
    res$species_name[3], paste(test$subfamily_name[3], "indet. species")
  )
  expect_true(is.na(res$family_name[2]))
  expect_true(is.na(res$subfamily_name[2]))
  expect_true(is.na(res$order_name[3]))

})
