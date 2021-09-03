#    test <- read.csv("tests/testthat/data_toclean_bold.csv")

test_that("Test refdb_clean_seq_remove_gaps", {

  test <- read.csv("data_toclean_bold.csv")
  test <- refdb_set_fields_BOLD(test)

  expect_equal(test$nucleotides[1], bioseq::dna("TAC---------TTT----TT"))
  res <- refdb_clean_seq_remove_gaps(test)
  expect_equal(res$nucleotides[1], bioseq::dna("TACTTTTT"))

})



test_that("Test refdb_clean_seq_remove_sideN", {

  test <- read.csv("data_toclean_bold.csv")
  test <- refdb_set_fields_BOLD(test)

  expect_equal(test$nucleotides[2], bioseq::dna("NNNNNNTACTTANNNNNN"))
  res <- refdb_clean_seq_remove_sideN(test, side = "left")
  expect_equal(res$nucleotides[2], bioseq::dna("TACTTANNNNNN"))
  res <- refdb_clean_seq_remove_sideN(test, side = "right")
  expect_equal(res$nucleotides[2], bioseq::dna("NNNNNNTACTTA"))
  res <- refdb_clean_seq_remove_sideN(test, side = "both")
  expect_equal(res$nucleotides[2], bioseq::dna("TACTTA"))

})



test_that("Test refdb_clean_seq_crop_primers", {

  test <- read.csv("data_toclean_bold.csv")
  test <- refdb_set_fields_BOLD(test)

  expect_equal(test$nucleotides[3], bioseq::dna("TGGTCGGCACTTCACTTAGTCTTTTTCTTATACCA"))
  res <- refdb_clean_seq_crop_primers(test,
                                      primer_forward = "GGCACT",
                                      primer_reverse = "TTTTT",
                                      max_error_in = 0,
                                      max_error_out = 0,
                                      include_primers = TRUE)
  expect_equal(res$nucleotides[3], bioseq::dna("GGCACTTCACTTAGTCTTTTT"))

})


test_that("Test refdb_clean_tax_remove_blank", {

  test <- read.csv("data_toclean_bold.csv")
  test <- refdb_set_fields_BOLD(test)

  expect_equal(test$subfamily_name[2], "Baetinae ")
  expect_equal(test$species_name[3], "Baetis  alpinus")
  res <- refdb_clean_tax_remove_blank(test)
  expect_equal(res$subfamily_name[2], "Baetinae")
  expect_equal(res$species_name[3], "Baetis alpinus")

})


test_that("Test refdb_clean_tax_remove_extra", {

  test <- read.csv("data_toclean_bold.csv")
  test <- refdb_set_fields_BOLD(test)

  test$species_name[1] <- "Cloeon dipterum ag415"
  res <- refdb_clean_tax_remove_extra(test)
  expect_equal(res$species_name[1], "Cloeon dipterum")

  test$species_name[1] <- "Cloeon dipterum CBG"
  res <- refdb_clean_tax_remove_extra(test)
  expect_equal(res$species_name[1], "Cloeon dipterum")

  test$species_name[1] <- "Cloeon dipterum X"
  res <- refdb_clean_tax_remove_extra(test)
  expect_equal(res$species_name[1], "Cloeon dipterum")

  test$species_name[1] <- "Cloeon dipterum x"
  res <- refdb_clean_tax_remove_extra(test)
  expect_equal(res$species_name[1], "Cloeon dipterum")

  test$species_name[1] <- "Cloeon dipterum (species)"
  res <- refdb_clean_tax_remove_extra(test)
  expect_equal(res$species_name[1], "Cloeon dipterum")

  test$species_name[1] <- "Cloeon (gen) dipterum"
  res <- refdb_clean_tax_remove_extra(test)
  expect_equal(res$species_name[1], "Cloeon dipterum")
})



test_that("Test refdb_clean_tax_harmonize_nomenclature", {

  test <- read.csv("data_toclean_bold.csv")
  test <- refdb_set_fields_BOLD(test)

  test$species_name[1] <- "Cloeon species"
  res <- refdb_clean_tax_harmonize_nomenclature(test)
  expect_equal(res$species_name[1], "Cloeon sp.")

  test$genus_name[1] <- "Baetidae familia"
  res <- refdb_clean_tax_harmonize_nomenclature(test)
  expect_equal(res$genus_name[1], "Baetidae fam.")

  test$species_name[1] <- "Cloeon confer dipterum"
  res <- refdb_clean_tax_harmonize_nomenclature(test)
  expect_equal(res$species_name[1], "Cloeon cf. dipterum")

})



test_that("Test refdb_clean_tax_remove_uncertainty", {

  test <- read.csv("data_toclean_bold.csv")
  test <- refdb_set_fields_BOLD(test)

  test$species_name[1] <- "Cloeon cf. dipterum"
  res <- refdb_clean_tax_remove_uncertainty(test)
  expect_equal(res$species_name[1], "Cloeon dipterum")

  test$species_name[1] <- "Cloeon aff. dipterum"
  res <- refdb_clean_tax_remove_uncertainty(test)
  expect_equal(res$species_name[1], "Cloeon dipterum")

  test$species_name[1] <- "Cloeon sp. inc. dipterum"
  res <- refdb_clean_tax_remove_uncertainty(test)
  expect_equal(res$species_name[1], "Cloeon dipterum")

})


test_that("Test refdb_clean_tax_remove_subsp", {

  test <- read.csv("data_toclean_bold.csv")
  test <- refdb_set_fields_BOLD(test)

  test$species_name[1] <- "Cloeon dipterum var. dipterum"
  res <- refdb_clean_tax_remove_subsp(test)
  expect_equal(res$species_name[1], "Cloeon dipterum")

})



test_that("Test refdb_clean_tax_NA", {

  test <- read.csv("data_toclean_bold.csv")
  test <- refdb_set_fields_BOLD(test)

  test$species_name[1] <- ""
  res <- refdb_clean_tax_NA(test)
  expect_equal(res$species_name[1], as.character(NA))

})
