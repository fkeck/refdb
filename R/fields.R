
# Use NA to ignore, use NULL to remove

refdb_set_fields <- function(x,
                             source = NA,
                             id = NA,
                             organism = NA,
                             taxonomy = NA,
                             sequence = NA,
                             marker = NA) {

  fields_list <- attr(x, "refdb_fields")
  if(is.null(fields_list)) {
    fields_list <- list()
  }

  error_field_col <- function(x, field, col) {
    if(!source %in% colnames(x)) {
      stop("The value ", col, " used for ", field,
           " does not match with any column of x.")
    }
  }


  if (!is.na(source)) {
    error_field_col(x, "source", source)
    fields_list$source <- source
    x[, source] <- as.character(x[, source, drop = TRUE])
  }

  if (!is.na(id)) {
    error_field_col(x, "id", id)
    fields_list$id <- id
    x[, id] <- as.character(x[, id, drop = TRUE])
  }

  if (!is.na(organism)) {
    error_field_col(x, "organism", organism)
    fields_list$organism <- organism
    x[, organism] <- as.character(x[, organism, drop = TRUE])
  }

  if (!all(is.na(taxonomy))) {

    if(!all(taxonomy %in% colnames(x))) {
      stop("The values: ", setdiff(taxonomy, colnames(x)),
           " used for taxonomy do not match with any column of x.")
    }

    if(!all(names(taxonomy) %in% valid_taxo_rank())) {
      stop("Taxonomic ranks are not valid. The `taxonomy` object must be a named vector. Use valid_taxo_rank() to find valid ranks.")
    }
    taxo_order <- match(names(taxonomy), valid_taxo_rank())
    fields_list$taxonomy <- taxonomy[order(taxo_order, na.last = NA)]
    x[, taxonomy] <- apply(x[, taxonomy], 2, as.character)
  }

  if (!is.na(sequence)) {
    error_field_col(x, "sequence", sequence)
    fields_list$sequence <- sequence
    x[, sequence] <- as.character(x[, sequence, drop = TRUE])
  }

  if (!is.na(marker)) {
    error_field_col(x, "marker", marker)
    fields_list$marker <- marker
    x[, marker] <- as.character(x[, marker, drop = TRUE])
  }

  attr(x, "refdb_fields") <- fields_list

  return(x)
}


# refdb_add_field <- function(x,
#                             field_name,
#                             col_name,
#                             col_value = NULL) {
#   if(is.null(col_value)){
#     stopifnot(col_name %in% colnames(x))
#   } else {
#     x[[col_name]] <- col_value
#   }
#
#   fields <- attr(x, "refdb_fields")
#   fields[[field_name]] <- col_name
#   attr(x, "refdb_fields") <- fields
#
#   return(x)
# }



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
         " -- See function `refdb_set_fields`.")
  }

  diff <- setdiff(unlist(fields[what]), colnames(x))
  if(length(diff) > 0) {
    stop("Some declared fields associated to non-existing columns: ", diff,
         " -- See function `refdb_set_fields`.")
  }

  invisible()
}
