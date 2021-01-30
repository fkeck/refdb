

#' Remove illegal characters from genetic sequences
#'
#' @param x a reference database with a defined sequence field.
#'
#' @return
#' A reference database.
#' @export
#'
refdb_clean_seq_remove_illegal <- function(x) {
  check_fields(x, "sequence")
  col <- attributes(x)$refdb$sequence
  x[[col]] <- bioseq::as_dna(x[[col]])
  return(x)
}


#' Remove gaps from genetic sequences
#'
#' @param x a reference database with a defined sequence field.
#'
#' @return
#' A reference database.
#' @export
#'
refdb_clean_seq_remove_gaps <- function(x) {
  check_fields(x, "sequence")
  col <- attributes(x)$refdb$sequence
  x[[col]] <- bioseq::seq_remove_pattern(x[[col]], "-")
  return(x)
}


#' Remove repeated side N from genetic sequences
#'
#' @param x a reference database with a defined sequence field.
#'
#' @return
#' A reference database.
#' @export
#'
refdb_clean_seq_remove_sideN <- function(x, side = "both") {

  check_fields(x, "sequence")
  col <- attributes(x)$refdb$sequence

  if(side == "in" | side == "both") {
    x[[col]] <- bioseq::seq_remove_pattern(x[[col]], "^N+")
  }

  if(side == "out" | side == "both") {
    x[[col]] <- bioseq::seq_remove_pattern(x[[col]], "N+$")
  }
  return(x)
}


#' Crop genetic sequences with a set of primers
#'
#' @param x a reference database with a defined sequence field.
#' @param pattern_in
#' @param pattern_out
#'
#' @return
#' A reference database.
#' @export
#'
refdb_clean_seq_crop_region <- function(x, pattern_in, pattern_out) {
  check_fields(x, "sequence")
  col <- attributes(x)$refdb$sequence

  x[[col]] <- bioseq::seq_crop_pattern(x[[col]],
                                       pattern_in = pattern_in,
                                       pattern_out = pattern_out)
  return(x)
}


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


#' Remove blank characters from taxonomic names
#'
#' @param x a reference database.
#' @param cols an optional vector of column names.
#' If \code{NULL} (default), the function is applied to the columns
#' associated with the \code{taxonomy} and \code{organism} fields.
#'
#' @return
#' A reference database.
#' @export
#'
refdb_clean_tax_remove_blank <- function(x, cols = NULL) {

  if(is.null(cols)) {
    cols <- c(attributes(x)$refdb$taxonomy,
             attributes(x)$refdb$organism)
  } else {
    # Maybe check for columns
  }

  x[, cols] <- apply(x[, cols], 2, stringr::str_trim)
  x[, cols] <- apply(x[, cols], 2, stringr::str_squish)
  return(x)
}


#' Remove extra words from taxonomic names
#'
#' @param x a reference database.
#' @param cols an optional vector of column names.
#' If \code{NULL} (default), the function is applied to the columns
#' associated with the \code{taxonomy} and \code{organism} fields.
#'
#' @return
#' A reference database.
#' @export
#'
refdb_clean_tax_remove_extra <- function(x, cols = NULL) {

  if(is.null(cols)) {
    cols <- c(attributes(x)$refdb$taxonomy,
              attributes(x)$refdb$organism)
  } else {
    # Maybe check for columns
  }

  .remove_fun <- function(x) {
    x <- stringr::str_remove_all(x, #Words containing numbers
                       "(?=\\S*['-])([a-zA-Z'-]+)\\d*(?=\\S*['-])([a-zA-Z'-]+)")
    x <- stringr::str_remove_all(x,
                       "\\w*[0-9]+\\w*\\s*")
    x <- stringr::str_remove_all(x, # Words in uppercase
                       "(?=\\S*['-])([a-zA-Z'-]+)[A-Z]{2,}(?=\\S*['-])([a-zA-Z'-]+)")
    x <- stringr::str_remove_all(x,
                       "\\w*[A-Z]{2,}\\w*\\s*")
    x <- stringr::str_remove_all(x, # Words of one character
                       "\\b\\w\\b")
    x <- stringr::str_remove_all(x, # Words between parentheses
                       "\\(.*\\)")
    return(x)
  }

  x[, cols] <- apply(x[, cols], 2, .remove_fun)
  return(x)
}


#' Harmonize taxonomic name nomenclature
#'
#' @param x a reference database.
#' @param cols an optional vector of column names.
#' If \code{NULL} (default), the function is applied to the columns
#' associated with the \code{taxonomy} and \code{organism} fields.
#'
#' @return
#' A reference database.
#' @export
#'
refdb_clean_tax_harmonize_nomenclature <- function(x, cols = NULL) {

  if(is.null(cols)) {
    cols <- c(attributes(x)$refdb$taxonomy,
              attributes(x)$refdb$organism)
  } else {
    # Maybe check for columns
  }

  .replace_fun <- function(x) {
    x <- stringr::str_replace(x, " familia ", " fam. ")
    x <- stringr::str_replace(x, " genus ", " gen. ")
    x <- stringr::str_replace(x, " species ", " sp. ")
    x <- stringr::str_replace(x, " subspecies | subsp. ", " ssp. ")

    x <- stringr::str_replace(x, " sp. plurimae ", " sp. pl. ")
    x <- stringr::str_replace(x, " g. sp. ", " gen. sp. ")

    x <- stringr::str_replace(x, " sp. indeterminabilis | sp. indeterminata | ind. | indet. | sp. ind. ", " sp. indet. ")
    x <- stringr::str_replace(x, " sp. nova | nova sp. | spec. nov. | sp. n. | nov. sp. | nov. spec. | n. sp. ", " sp. nov. ")

    x <- stringr::str_replace(x, " sp. affinis | sp. aff. ", " aff. ")
    x <- stringr::str_replace(x, " sp. proxima | sp. prox. ", " prox. ")
    x <- stringr::str_replace(x, " sp. near | sp. nr. ", " nr. ")
    x <- stringr::str_replace(x, " ex grege | gr. ", " ex gr. ")
    x <- stringr::str_replace(x, " confer | cfr. | conf. | sp. cf. ", " cf. ")
    x <- stringr::str_replace(x, " sp. incerta | inc. | \\? ", " sp. inc. ")
    x <- stringr::str_replace(x, " stetit ", " stet. ")

    x <- stringr::str_replace(x, " sensu lato ", " s. l. ")
    x <- stringr::str_replace(x, " incertae sedis ", " inc. sed. ")

    return(x)
  }

  x[, cols] <- apply(x[, cols], 2, .replace_fun)
  return(x)

}


#' Remove terms indicating uncertainty in taxonomic names
#'
#' @param x a reference database.
#' @param cols an optional vector of column names.
#' If \code{NULL} (default), the function is applied to the columns
#' associated with the \code{taxonomy} and \code{organism} fields.
#'
#' @return
#' A reference database.
#' @export
#'
refdb_clean_tax_remove_uncertainty <- function(x, cols = NULL) {

  if(is.null(cols)) {
    cols <- c(attributes(x)$refdb$taxonomy,
              attributes(x)$refdb$organism)
  } else {
    # Maybe check for columns
  }

  x[, cols] <- apply(x[, cols], 2,
                     stringr::str_replace,
                     pattern = " aff\\. | cf\\. | prox\\. | nr\\. | sp\\. inc\\. ",
                     replacement = " ")

  return(x)
}


#' Remove subspecific information from taxonomic names
#'
#' @param x a reference database.
#' @param cols an optional vector of column names.
#' If \code{NULL} (default), the function is applied to the columns
#' associated with the \code{taxonomy} and \code{organism} fields.
#'
#' @return
#' A reference database.
#' @export
#'
refdb_clean_tax_remove_subsp <- function(x, cols = NULL) {

  if(is.null(cols)) {
    cols <- c(attributes(x)$refdb$taxonomy,
              attributes(x)$refdb$organism)
  } else {
    # Maybe check for columns
  }

  x[, cols] <- apply(x[, cols], 2,
                     stringr::str_replace,
                     pattern = " var\\. .*$| v\\. .*$| varietas .*$| forma .*$| f\\. .*$| morph .*$| form .*$| biotype .*$| isolate .*$| pathogroup .*$| serogroup .*$| serotype .*$| strain .*$| aberration .*$| abberatio .*$| ab\\. .*$",
                     replacement = "")

  return(x)

}

#' Convert missing taxonomic names to NA
#'
#' @param x a reference database.
#' @param cols an optional vector of column names.
#' If \code{NULL} (default), the function is applied to the columns
#' associated with the \code{taxonomy} and \code{organism} fields.
#'
#' @return
#' A reference database.
#' @export
#'
refdb_clean_tax_NA <- function(x, cols = NULL) {

  if(is.null(cols)) {
    cols <- c(attributes(x)$refdb$taxonomy,
              attributes(x)$refdb$organism)
  } else {
    # Maybe check for columns
  }

  .replace_fun <- function(x) {
    x <- ifelse(x == "", NA, x)
    x <- ifelse(stringr::str_detect(x, " stetit$| stet\\.| sp[0-9]| sp\\.| sp$| sp "), NA, x)
  }

  x[, cols] <- apply(x[, cols], 2, .replace_fun)
  return(x)

}
