


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
#' @param side which side to clean.
#' Can be one of \code{"left"}, \code{"right"} or \code{"both"} (default).
#'
#' @return
#' A reference database.
#' @export
#'
refdb_clean_seq_remove_sideN <- function(x, side = "both") {

  check_fields(x, "sequence")
  col <- attributes(x)$refdb$sequence

  if(side == "left" | side == "both") {
    x[[col]] <- bioseq::seq_remove_pattern(x[[col]], "^N+")
  }

  if(side == "right" | side == "both") {
    x[[col]] <- bioseq::seq_remove_pattern(x[[col]], "N+$")
  }
  return(x)
}


#' Crop genetic sequences with a set of primers
#'
#' @param x a reference database with a defined sequence field.
#' @param primer_forward primer forward.
#' @param primer_reverse primer reverse.
#' @param max_error_in,max_error_out maximum error for a match
#' (frequency based on primer length).
#' @param include_primers a logical indicating whether the detected primers are
#' included in the cropped sequences.
#'
#' @return
#' A reference database.
#' @export
#'
refdb_clean_seq_crop_primers <- function(x,
                                         primer_forward,
                                         primer_reverse,
                                         max_error_in = 0.1,
                                         max_error_out = 0.1,
                                         include_primers = TRUE) {

  check_fields(x, "sequence")
  col <- attributes(x)$refdb$sequence

  x[[col]] <- bioseq::seq_crop_pattern(x[[col]],
                                       pattern_in = primer_forward,
                                       pattern_out = primer_reverse,
                                       max_error_in = max_error_in,
                                       max_error_out = max_error_out,
                                       include_patterns = include_primers)
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

  x <- sanity_check(x, cols = cols, clean_whitespace = FALSE)

  return(x)
}


#' Remove extra words from taxonomic names
#'
#' @param x a reference database.
#' @param cols an optional vector of column names.
#' If \code{NULL} (default), the function is applied to the columns
#' associated with the \code{taxonomy} and \code{organism} fields.
#' @details
#' As the function can match words like "g.", "s." or "x", which can
#' have a signification in some nomenclatures, it is recommended to
#' execute \link{refdb_clean_tax_harmonize_nomenclature} first.
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
    x <- stringr::str_remove_all(x, # Remove remaining dots
                       "\\.{2,}$| \\.+||^\\.+")
    x <- stringr::str_remove_all(x, # Remove word starting with uppercase (and remain string) after a word in lowercase/dot
                                 "(?<= [a-z\\.]{1,30}) [A-Z].*")
    x <- stringr::str_remove_all(x, # Remove words after sp. but not if the first is ending with a dot
                                 "(?<= sp\\.) [a-z]+($| ).*")
    x <- stringr::str_remove_all(x, # Remove word after sp. nov.
                                 "(?<=[A-Z][a-z]{2,30} [a-z]{2,30} sp\\. nov\\.).*")

    x <- stringr::str_squish(x)
    x <- stringr::str_trim(x)
    return(x)
  }

  x[, cols] <- apply(x[, cols], 2, .remove_fun)

  x <- sanity_check(x, cols = cols)

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
    x <- stringr::str_replace(x, " familia($|\\s)", " fam. ")
    x <- stringr::str_replace(x, " genus($|\\s)", " gen. ")
    x <- stringr::str_replace(x, " species($|\\s)", " sp. ")
    x <- stringr::str_replace(x, " subspecies($|\\s)| subsp\\.($|\\s)", " ssp. ")

    x <- stringr::str_replace(x, " sp\\. plurimae($|\\s)", " sp. pl. ")
    x <- stringr::str_replace(x, " g\\. sp\\.($|\\s)", " gen. sp. ")
    x <- stringr::str_replace(x, " sp\\. indeterminabilis($|\\s)| sp\\. indeterminata($|\\s)| ind\\.($|\\s)| indet\\.($|\\s)| sp\\. ind\\.($|\\s)", " sp. indet. ")
    x <- stringr::str_replace(x, " sp\\. nova($|\\s)| nova sp\\.($|\\s)| spec\\. nov\\.($|\\s)| sp\\. n\\.($|\\s)| nov\\. sp\\.($|\\s)| nov\\. spec\\.($|\\s)| n\\. sp\\.($|\\s)", " sp. nov. ")

    x <- stringr::str_replace(x, " sp\\. affinis($|\\s)| sp\\. aff\\.($|\\s)", " aff. ")
    x <- stringr::str_replace(x, " sp\\. proxima($|\\s)| sp\\. prox\\.($|\\s)", " prox. ")
    x <- stringr::str_replace(x, " sp\\. near($|\\s)| sp\\. nr\\.($|\\s)", " nr. ")
    x <- stringr::str_replace(x, " confer($|\\s)| cfr\\.($|\\s)| conf\\.($|\\s)| sp\\. cf\\.($|\\s)", " cf. ")
    x <- stringr::str_replace(x, " sp\\. incerta($|\\s)| inc\\.($|\\s)| \\?($|\\s)", " sp. inc. ")
    x <- stringr::str_replace(x, " stetit($|\\s)", " stet. ")

    x <- stringr::str_replace(x, " sensu lato($|\\s)", " sen. lat. ")
    x <- stringr::str_replace(x, " incertae sedis($|\\s)", " inc. sed. ")

    x <- stringr::str_replace(x, "(^|\\s)sp\\. complex($|\\s)|(^|\\s)complex sp\\.($|\\s)|(^|\\s)group($|\\s)|(^|\\s)group sp\\.($|\\s)|(^|\\s)aggregate sp\\.($|\\s)|(^|\\s)ex grege($|\\s)|(^|\\s)gr\\.($|\\s)|(^|\\s)ex gr\\.($|\\s)", " complex ")
    x <- stringr::str_replace(x, "(^|\\s)complex sp\\.($|\\s)", " complex ")

    x <- stringr::str_replace_all(x, .REGEX_HYBRID, " * ")

    return(x)
  }

  x[, cols] <- apply(x[, cols], 2, .replace_fun)

  x <- sanity_check(x, cols = cols)

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
                     pattern = .REGEX_UNCERTAIN,
                     replacement = " ")

  x <- sanity_check(x, cols = cols)

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
                     pattern = .REGEX_SUBSP,
                     replacement = "")

  x <- sanity_check(x, cols = cols)

  return(x)

}

# TODO
# Remove taxonomic qualifiers eg. sp. nov.
# and everything following
# To be run after harmonizing nomenclature (or maybe enforce it inside the function)
# refdb_clean_tax_remove_qualifiers <- function(x, cols = NULL) {
#
#   if(is.null(cols)) {
#     cols <- c(attributes(x)$refdb$taxonomy,
#               attributes(x)$refdb$organism)
#   } else {
#     # Maybe check for columns
#   }
#
#   x[, cols] <- apply(x[, cols], 2,
#                      stringr::str_replace,
#                      pattern = "",
#                      replacement = "")
#
#   x <- sanity_check(x, cols = cols)
#
#   return(x)
#
# }

#' Convert missing taxonomic names to NA
#'
#' @param x a reference database.
#' @param cols an optional vector of column names.
#' If \code{NULL} (default), the function is applied to the columns
#' associated with the \code{taxonomy} and \code{organism} fields.
#' @param hybrid hybrids are converted to NA (default \code{TRUE}).
#' @param uncertain taxa with qualifiers of uncertainty (cf., aff., etc.)
#' are converted to NA (default \code{FALSE}).
#'
#' @return
#' A reference database.
#' @export
#'
refdb_clean_tax_NA <- function(x, cols = NULL, hybrid = TRUE, uncertain = FALSE) {

  if(is.null(cols)) {
    cols <- c(attributes(x)$refdb$taxonomy,
              attributes(x)$refdb$organism)
  } else {
    # Maybe check for columns
  }

  rgx <- .REGEX_NOT_ID

  if(hybrid) {
    rgx <- paste0(rgx, "|", .REGEX_HYBRID)
  }

  if(uncertain) {
    rgx <- paste0(rgx, "|", .REGEX_UNCERTAIN)
  }

  .replace_fun <- function(x) {
    x <- ifelse(x == "", NA, x)
    x <- ifelse(stringr::str_detect(x, rgx), NA, x)
  }

  x[, cols] <- apply(x[, cols], 2, .replace_fun)

  if("species" %in% names(cols)) {
    x[[cols["species"]]][!str_detect(x[[cols["species"]]], "[A-Za-z] +[A-Za-z]")] <- NA
  }

  return(x)
}




# refdb_clean_tax_smart <- function(x, cols = NULL) {
#   refdb_clean_tax_harmonize_nomenclature
#   refdb_clean_tax_remove_extra
#   refdb_clean_tax_NA
# }




sanity_check <- function(x, cols,
                         clean_whitespace = TRUE,
                         clean_na = TRUE) {

  if(clean_whitespace) {
    x[, cols] <- apply(x[, cols], 2, stringr::str_trim)
    x[, cols] <- apply(x[, cols], 2, stringr::str_squish)
  }

  if(clean_na) {
    x[, cols] <- apply(x[, cols], 2, function(x) ifelse(x == "", NA, x))
  }

  return(x)
}



#### REGEXES ####

.REGEX_UNCERTAIN <- " aff\\. | cf\\. | prox\\. | nr\\. | sp\\. inc\\. |^aff\\.|^cf\\.|^prox\\.|^nr\\.|^sp\\. inc\\.| aff\\.$| cf\\.$| prox\\.$| nr\\.$| sp\\. inc\\.$"
.REGEX_NOT_ID <- "(^|\\s)stetit($|\\s)|(^|\\s)stet\\.($|\\s)|(^|\\s)sp[0-9]($|\\s)|([A-Z][a-z]+ sp\\.)| gen\\. sp\\.|(^|\\s)unclassified($|\\s)|^[:blank:]+$"
.REGEX_HYBRID <- " \\* | \\u00D7 | [xX] (?=[A-Z])"
.REGEX_SUBSP <- " var\\. .*$| v\\. .*$| varietas .*$| forma .*$| f\\. .*$| morph .*$| form .*$| biotype .*$| isolate .*$| pathogroup .*$| serogroup .*$| serotype .*$| strain .*$| aberration .*$| abberatio .*$| ab\\. .*$"


