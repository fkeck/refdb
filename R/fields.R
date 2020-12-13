

refdb_set_fields <- function(x,
                             source = NULL,
                             id = NULL,
                             organism = NULL,
                             taxonomy = NULL,
                             sequence = NULL,
                             marker = NULL) {

  stopifnot(source %in% colnames(x))
  stopifnot(id %in% colnames(x))
  stopifnot(organism %in% colnames(x))
  stopifnot(all(taxonomy %in% colnames(x)))
  stopifnot(sequence %in% colnames(x))
  stopifnot(marker %in% colnames(x))

  attr(x, "refdb_fields") <- list(source = source,
                                  id = id,
                                  organism = organism,
                                  taxonomy = taxonomy,
                                  sequence = sequence,
                                  marker = marker)
  return(x)
}


refdb_add_field <- function(x,
                            field_name,
                            col_name,
                            col_value = NULL) {
  if(is.null(col_value)){
    stopifnot(col_name %in% colnames(x))
  } else {
    x[[col_name]] <- col_value
  }

  fields <- attr(x, "refdb_fields")
  fields[[field_name]] <- col_name
  attr(x, "refdb_fields") <- fields

  return(x)
}



#' Internal check for fields
#'
#' @param x a reference database (tibble object).
#' @param what a vector of fields to be checked.
#'
#' @return
#' Invisible or error.
#'
check_fields <- function(x,
                         what = c("source", "id", "organism",
                                  "taxonomy", "sequence", "marker")) {
  fields <- attributes(x)$refdb_fields

  diff <- setdiff(names(fields), what)
  if(length(diff) > 0) {
    stop("Missing field: ", diff,
         " -- See function `refdb_set_fields` and `refdb_add_field`.")
  }

  diff <- setdiff(unlist(fields[what]), colnames(x))
  if(length(diff) > 0) {
    stop("Some declared fields associated to non-existing columns: ", diff,
         " -- See function `refdb_set_fields` and `refdb_add_field`.")
  }

}
