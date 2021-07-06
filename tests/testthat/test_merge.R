
# test_2 <- read.csv("tests/testthat/data_fragilaria_ncbi.csv")
# test_1 <- read.csv("tests/testthat/data_baetidae_bold.csv")



test_that("Test merging databases", {
  test_1 <- read.csv("data_baetidae_bold.csv")
  test_1 <- refdb_set_fields_BOLD(test_1)

  test_2 <- read.csv("data_fragilaria_ncbi.csv")
  test_2 <- refdb_set_fields_NCBI(test_2)

  res <- refdb_merge(test_1, test_2, keep = "fields_all")
  expect_equal(nrow(res), nrow(test_1) + nrow(test_2))
  expect_equal(ncol(res),
               length(unique(c(names(unlist(attributes(test_1)$refdb_fields)),
                               names(unlist(attributes(test_2)$refdb_fields))))))
})


test_that("Test merging databases", {
  test_1 <- read.csv("data_baetidae_bold.csv")
  test_1 <- refdb_set_fields_BOLD(test_1)

  test_2 <- read.csv("data_fragilaria_ncbi.csv")
  test_2 <- refdb_set_fields_NCBI(test_2)

  res <- refdb_merge(test_1, test_2, keep = "fields_shared")
  expect_equal(nrow(res), nrow(test_1) + nrow(test_2))
  expect_equal(ncol(res),
               length(intersect(names(unlist(attributes(test_1)$refdb_fields)),
                               names(unlist(attributes(test_2)$refdb_fields)))))
})


test_that("Test merging databases", {
  test_1 <- read.csv("data_baetidae_bold.csv")
  test_1 <- refdb_set_fields_BOLD(test_1)

  test_2 <- read.csv("data_fragilaria_ncbi.csv")
  test_2 <- refdb_set_fields_NCBI(test_2)

  res <- refdb_merge(test_1, test_2, keep = "all")
  expect_equal(nrow(res), nrow(test_1) + nrow(test_2))
  expect_equal(ncol(res), 26L)
})
