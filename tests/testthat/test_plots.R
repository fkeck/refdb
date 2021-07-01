# test <- read.csv("tests/testthat/data_baetidae_bold.csv")

test_that("Test plot functions", {
  test <- read.csv("data_baetidae_bold.csv")
  test <- refdb_set_fields_BOLD(test)

  res <- igraph_from_taxo(test)

  expect_s3_class(res, class = "igraph")
  expect_length(res, 10L)

  p <- refdb_plot_tax_tree(test)
  p <- refdb_plot_tax_treemap(test)
  p <- refdb_plot_seqlen_hist(test)

})
