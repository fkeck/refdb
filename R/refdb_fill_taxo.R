


#' Fill missing data in taxonomy
#'
#' Replace NA values in taxonomic classification
#' using downstream ranks.
#'
#' @param x a reference database.
#' @param qualifier a string to add the new labels.
#' Default ensure that \code{refdb_clean_tax_NA} will correctly
#' identify the label as NA.
#'
#' @return
#' A reference database.
#'
#' @seealso refdb_fill_tax_downstream to replace terminal NA values
#' using upstream data.
#'
#' @export
#'
refdb_fill_tax_upstream <- function(x, qualifier = "undef.") {
  check_fields(x, what = "taxonomy")
  col_tax <- attributes(x)$refdb_fields$taxonomy
  tax_mat <- x[, col_tax]
  tax_levels <- names(col_tax)

  fun_row <- function(row) {
    row_len <- length(row)
    for(i in seq(row_len - 1)) {
      if(is.na(row[i])) {
        repl <- stats::na.exclude(row[seq(i + 1, row_len)])
        if(length(repl) == 0) {
          row[i] <- NA
        } else {
          row[i] <- paste(repl[1], qualifier, tax_levels[i])
        }
      }
    }
    return(row)
  }

  res <- apply(tax_mat, 1, function(y) {
    if(any(is.na(y))) {
      y <- fun_row(y)
    }
    return(y)
  })

  x[, col_tax] <- t(res)
  return(x)
}


#' Fill missing data in taxonomy
#'
#' Replace NA values in taxonomic classification
#' using upstream ranks.
#'
#' @param x a reference database.
#' @param qualifier a string to add the new labels.
#' Default ensure that \code{refdb_clean_tax_NA} will correctly
#' identify the label as NA.
#'
#' @return
#' A reference database.
#'
#' @seealso refdb_fill_tax_upstream to replace NA values
#' using downstream data.
#'
#' @export
#'
refdb_fill_tax_downstream <- function(x, qualifier = "indet.") {
  check_fields(x, what = "taxonomy")
  col_tax <- attributes(x)$refdb_fields$taxonomy
  tax_mat <- x[, col_tax]
  tax_levels <- names(col_tax)

  fun_row <- function(row) {
    row_len <- length(row)
    idx_max_not_na <- max(which(!is.na(row)))
    if(idx_max_not_na == row_len) {
      return(row)
    } else {
      idx_na <- seq(idx_max_not_na + 1, row_len)
      row[idx_na] <- paste(row[idx_max_not_na], qualifier, tax_levels[idx_na])
      return(row)
    }
  }

  res <- apply(tax_mat, 1, function(y) {
    if(any(is.na(y))) {
      y <- fun_row(y)
    }
    return(y)
  })

  x[, col_tax] <- t(res)
  return(x)
}

