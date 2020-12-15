
#' Merge reference databases
#'
#' Merge several reference database by common fields.
#'
#' @param ... reference databases (tibbles).
#' @param keep determines which columns to keep.
#' Can be \code{"fields_all"} (default), \code{"fields_shared"}
#' or \code{"all"} (see Details).
#'
#' @details
#'
#' Columns are merged only if they are associated to the same field.
#'
#' The \code{keep} argument determines which columns are returned as follow.
#' \code{"fields_all"} (the default) returns all the fields
#' existing in all the reference databases.
#' \code{"fields_shared"} returns only the fields shared by
#' all the reference databases.
#' \code{"all"} returns all the columns of all the databases.
#' Columns which are not associated to a field are not merged and are prefixed
#' with the name of the object they originated from.
#'
#'
#' @return a merged reference database (tibble).
#' @export
#'
refdb_merge <-function(..., keep = "fields_all") {

  x <- list(...)

  stopifnot("All reference db objects must be dataframes" =
              all(sapply(x, is.data.frame))
            )
  lapply(x, check_fields)


  fields_list <- lapply(x, function(x) unlist(attributes(x)$refdb_fields))
  fields_list_names <- lapply(fields_list, names)

  fields_union <- Reduce(union, fields_list_names)
  fields_inter <- Reduce(intersect, fields_list_names)

  if(keep == "fields_shared") {
    fields_keep <- lapply(x, function(x) unlist(attributes(x)$refdb_fields)[fields_inter])

    out <- mapply(function(x, y) {
      res <- x[, y]
      colnames(res) <- fields_keep[[1]]
      return(res)
    },
    x = x,
    y = fields_keep,
    SIMPLIFY = FALSE)
  }


  if(keep == "fields_all") {
    fields_keep <- lapply(fields_list, function(x) x[intersect(names(x), fields_union)])

    out <- mapply(function(x, y) {
      res <- x[, y]
      colnames(res) <- unlist(fields_list)[names(y)]
      return(res)
    },
    x = x,
    y = fields_keep,
    SIMPLIFY = FALSE)
  }


  if(keep == "all") {
    fields_keep <- lapply(fields_list, function(x) x[intersect(names(x), fields_union)])

    out <- mapply(function(x, y) {
      res <- x[, c(y, setdiff(colnames(x), y))]
      colnames(res) <- c(unlist(fields_list)[names(y)],
                         paste0("AAA_", setdiff(colnames(x), y)))
      return(res)
    },
    x = x,
    y = fields_keep,
    SIMPLIFY = FALSE)
  }

  out <- dplyr::bind_rows(out)

  return(out)
}


