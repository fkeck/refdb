
context("Read external")

# test <- read.csv("tests/testthat/data_fragilaria_ncbi.csv")
# test <- read.csv("tests/testthat/data_baetidae_bold.csv")


test_that("Test manually setting fields", {
  test <- read.csv("data_fragilaria_ncbi.csv")

  test <- refdb_set_fields(test,
                   source = "source",
                   id = "id",
                   taxonomy = c(superkingdom = "superkingdom",
                                kingdom = "kingdom",
                                phylum = "phylum",
                                subphylum = "subphylum",
                                class = "class",
                                subclass = "subclass",
                                infraclass = "infraclass",
                                order = "order",
                                suborder = "suborder",
                                infraorder = "infraorder",
                                superfamily = "superfamily",
                                family = "family",
                                genus = "genus",
                                species = "species"),
                   sequence = "sequence",
                   marker = "gene")

  expect_s3_class(test, "data.frame")
  expect_s3_class(test$sequence, "bioseq_dna")
  expect_equal(ncol(test), 20L)
  expect_equal(nrow(test), 41L)
  expect_length(attributes(test)$refdb_fields, 5L)
  expect_equal(attributes(test)$refdb_fields,
               list(source = "source",
                    id = "id",
                    taxonomy = c(superkingdom = "superkingdom",
                                 kingdom = "kingdom", phylum = "phylum", subphylum = "subphylum",
                                 class = "class", subclass = "subclass", infraclass = "infraclass",
                                 order = "order", suborder = "suborder", infraorder = "infraorder",
                                 superfamily = "superfamily", family = "family", genus = "genus",
                                 species = "species"),
                    sequence = "sequence",
                    marker = "gene"))

})


test_that("Test setting fields with YAML", {
  test <- read.csv("data_baetidae_bold.csv")

  test <- refdb_set_fields(test, config_yaml = "data_yaml_bold.yml")

  expect_s3_class(test, "data.frame")
  expect_s3_class(test$nucleotides, "bioseq_dna")
  expect_equal(ncol(test), 16L)
  expect_equal(nrow(test), 30L)
  expect_length(attributes(test)$refdb_fields, 5L)
  expect_equal(attributes(test)$refdb_fields,
               list(source = "source",
                    id = "sequenceID",
                    taxonomy = c(phylum = "phylum_name",
                                 class = "class_name", order = "order_name", family = "family_name",
                                 subfamily = "subfamily_name", genus = "genus_name", species = "species_name"),
                    sequence = "nucleotides",
                    marker = "markercode"))
})




test_that("Test error catching in setting fields", {
  test <- read.csv("data_fragilaria_ncbi.csv")

  expect_error(refdb_set_fields(test,
                           taxonomy = c(superkingdom = "superkingdom",
                                        subphylum = "subphylum",
                                        class = "does_not_exist",
                                        subclass = "subclass"
                                        )),
               regexp = "do not match with any column")

  expect_error(refdb_set_fields(test,
                                taxonomy = c(superkingdom = "superkingdom",
                                             kingdom = "kingdom",
                                             weird_rank = "phylum",
                                             subphylum = "subphylum"
                                )),
               regexp = "Taxonomic ranks are not valid")

  expect_error(refdb_set_fields(test,
                                reference = "doesnotexist"),
               regexp = "do not match with any column")
  expect_error(refdb_set_fields(test,
                                source = "doesnotexist"),
               regexp = "does not match with any column")
  refdb_set_fields(test, organism = "species")
})



test_that("Test setting fields with BOLD function", {
  test <- read.csv("data_baetidae_bold.csv")

  test <- refdb_set_fields_BOLD(test)

  expect_s3_class(test, "data.frame")
  expect_s3_class(test$nucleotides, "bioseq_dna")
  expect_equal(ncol(test), 16L)
  expect_equal(nrow(test), 30L)
  expect_length(attributes(test)$refdb_fields, 5L)
  expect_equal(attributes(test)$refdb_fields,
               list(source = "source",
                    id = "sequenceID",
                    taxonomy = c(phylum = "phylum_name",
                                 class = "class_name", order = "order_name", family = "family_name",
                                 subfamily = "subfamily_name", genus = "genus_name", species = "species_name"),
                    sequence = "nucleotides",
                    marker = "markercode"))
})


test_that("Account NCBI ranks", {
  res <- ncbi_taxo_rank()
  expect_length(res, 13L)
})
