#    test <- read.csv("tests/testthat/data_tofilter_bold.csv")

test_that("Test refdb_filter_seq_length", {

  test <- read.csv("data_tofilter_bold.csv")
  test <- refdb_set_fields_BOLD(test)

  test$nucleotides[1] <- bioseq::dna("ACGT")
  res <- refdb_filter_seq_length(test, min_len = 5)
  expect_equal(nrow(res), nrow(test) - 1)

  test$nucleotides[1] <- bioseq::dna(paste(rep("ACGT", 10000), collapse = ""))
  res <- refdb_filter_seq_length(test, max_len = 1000)
  expect_equal(nrow(res), nrow(test) - 1)

  test$nucleotides[1] <- bioseq::dna("ACGT")
  test$nucleotides[2] <- bioseq::dna(paste(rep("ACGT", 10000), collapse = ""))
  res <- refdb_filter_seq_length(test, min_len = 5, max_len = 1000)
  expect_equal(nrow(res), nrow(test) - 2)
})

test_that("Test refdb_filter_seq_ambiguous", {

  test <- read.csv("data_tofilter_bold.csv")
  test <- refdb_set_fields_BOLD(test)

  test$nucleotides[1] <- bioseq::dna("ACGTACNNGTACGT")
  res <- refdb_filter_seq_ambiguous(test, max_ambig = 3)
  expect_equal(nrow(res), nrow(test))

  test$nucleotides[1] <- bioseq::dna("ACGTACNNNNGTACGT")
  res <- refdb_filter_seq_ambiguous(test, max_ambig = 3)
  expect_equal(nrow(res), nrow(test) - 1)

  test$nucleotides[1] <- bioseq::dna("ACGTACNNGTACRRRGT")
  res <- refdb_filter_seq_ambiguous(test, max_ambig = 3, char = c("N", "R"))
  expect_equal(nrow(res), nrow(test) - 1)

})


test_that("Test refdb_filter_seq_homopolymers", {

  test <- read.csv("data_tofilter_bold.csv")
  test <- refdb_set_fields_BOLD(test)

  test$nucleotides[1] <- bioseq::dna("ACGTACAAAAAAAAAAAAAAAAGTACGT")
  res <- refdb_filter_seq_homopolymers(test, max_len = 14)
  expect_equal(nrow(res), nrow(test) - 1)

  test$nucleotides[1] <- bioseq::dna("ACGTACAAAAAAAAAAAAAAAAGTACGT")
  res <- refdb_filter_seq_homopolymers(test, max_len = 20)
  expect_equal(nrow(res), nrow(test))

})



test_that("Test refdb_filter_seq_duplicates", {

  test <- read.csv("data_tofilter_bold.csv")
  test <- refdb_set_fields_BOLD(test)

  test[2, ] <- test[10, ]
  res <- refdb_filter_seq_duplicates(test)
  expect_equal(nrow(res), nrow(test) - 1)

})


test_that("Test refdb_filter_seq_stopcodon", {

  test <- read.csv("data_tofilter_bold.csv")
  test <- refdb_set_fields_BOLD(test)

  test$nucleotides[1] <- bioseq::dna("AAA")
  res <- refdb_filter_seq_stopcodon(test, code = 5, codon_frame = 1)
  expect_equal(nrow(res), 1)

  test$nucleotides[1] <- bioseq::dna("AAA")
  res <- refdb_filter_seq_stopcodon(test, code = 5, codon_frame = NA)
  expect_equal(nrow(res), nrow(test))

})


test_that("Test refdb_filter_seq_primer", {

  test <- read.csv("data_tofilter_bold.csv")
  test <- refdb_set_fields_BOLD(test)

  res <- refdb_filter_seq_primer(test,
                                 primer_forward = NULL,
                                 primer_reverse = NULL)
  expect_equal(nrow(res), nrow(test))


  test$nucleotides[1] <- bioseq::dna("ACGTAAAAAAACGTAACGTGGGGGGGGTA")
  test$nucleotides[10] <- bioseq::dna("ACGTAAAAAAACGTA")
  res <- refdb_filter_seq_primer(test,
                                 primer_forward = "AAAAAAA",
                                 primer_reverse = "GGGGGGGG")
  expect_equal(nrow(res), 1)

  res <- refdb_filter_seq_primer(test,
                                 primer_reverse = "GGGGGGGG")
  expect_equal(nrow(res), 1)

  res <- refdb_filter_seq_primer(test,
                                 primer_forward = "TTTTTTTTTTT")
  expect_equal(nrow(res), 0)
})



test_that("Test refdb_filter_tax_precision", {

  test <- read.csv("data_tofilter_bold.csv")
  test <- refdb_set_fields_BOLD(test)

  res <- refdb_filter_tax_precision(test, min_tax = "species_name")
  expect_equal(nrow(res), nrow(test))
  res <- refdb_filter_tax_precision(test, min_tax = "genus_name")
  expect_equal(nrow(res), nrow(test))

  test$species_name[3] <- NA
  res <- refdb_filter_tax_precision(test, min_tax = "species_name")
  expect_equal(nrow(res), nrow(test) - 1)

  test$species_name[3] <- "Cloeon dipterum"
  test$genus_name[3] <- NA
  res <- refdb_filter_tax_precision(test, min_tax = "species_name")
  expect_equal(nrow(res), nrow(test))
})


test_that("Test refdb_filter_ref_scope", {

  test <- read.csv("data_tofilter_bold.csv")
  test <- refdb_set_fields_BOLD(test)
  test$refs <- rep("REF_1", nrow(test))
  test <- refdb_set_fields(test, reference = "refs")

  res <- refdb_filter_ref_scope(test, max_tax = "genus_name")
  expect_equal(nrow(res), 0)

  res <- refdb_filter_ref_scope(test, max_tax = "family_name")
  expect_equal(nrow(res), nrow(test))

  test$class_name[3] <- "Anotherclass"
  test[3, attributes(test)$refdb_fields$taxonomy] <-
    rep("Another", length(attributes(test)$refdb_fields$taxonomy))

  res <- refdb_filter_ref_scope(test, max_tax = "family_name")
  expect_equal(nrow(res), 0)

  test <- read.csv("data_tofilter_bold.csv")
  test <- refdb_set_fields_BOLD(test)
  test$refs <- rep("REF_1", nrow(test))
  test$refs[1:10] <- rep("REF_2", 10)
  test <- refdb_set_fields(test, reference = "refs")

  res <- refdb_filter_ref_scope(test, max_tax = "genus_name")
  expect_equal(nrow(res), nrow(test) - 10)
})



test_that("Test refdb_filter_tax_na", {

  test <- read.csv("data_tofilter_bold.csv")
  test <- refdb_set_fields_BOLD(test)
  test[3, c("genus_name", "species_name")] <- NA
  res <- refdb_filter_tax_na(test)
  expect_equal(nrow(res), nrow(test) - 1)

  test <- read.csv("data_tofilter_bold.csv")
  test <- refdb_set_fields_BOLD(test)
  test[5, c("genus_name", "species_name")] <- NA
  test[5, c("genus_name")] <- "Neogen"
  res <- refdb_filter_tax_na(test)
  expect_equal(nrow(res), nrow(test))

  test <- read.csv("data_tofilter_bold.csv")
  test <- refdb_set_fields_BOLD(test)
  test[7, c("order_name")] <- NA
  res <- refdb_filter_tax_na(test)
  expect_equal(nrow(res), nrow(test))

  test <- read.csv("data_tofilter_bold.csv")
  test <- refdb_set_fields_BOLD(test)
  test[1:2, c("species_name")] <- NA
  res <- refdb_filter_tax_na(test)
  expect_equal(nrow(res), nrow(test))
})

