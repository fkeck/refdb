

#' Sample records within taxa
#'
#' This function can be useful to keep a maximum of records per taxa.
#'
#' @param x a reference database.
#' @param cols an optional vector of column names.
#' If \code{NULL} (default), the function is applied to the columns
#' associated with the \code{taxonomy} field.
#' @param n_max maximum number of records to keep for each taxa.
#'
#' @return
#' A reference database.
#' @export
#'
refdb_sample_tax <- function(x, n_max = 10, cols = NULL){
  check_fields(x, "taxonomy")

  if(is.null(cols)) {
    cols <- attributes(x)$refdb$taxonomy
  } else {
    # Maybe check for columns
  }

  x <- dplyr::group_by(x, !!!rlang::syms(unname(cols)))
  x <- dplyr::slice_sample(x, n = n_max)
  return(x)
}





