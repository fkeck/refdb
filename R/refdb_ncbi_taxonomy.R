


#' Replace the current taxonomy using the NCBI Taxonomy database
#'
#' @param x a reference database (tibble) with one or several columns
#' giving the taxonomy of each record and explicitly
#' indicated in the field taxonomy. See \link{refdb_set_fields}.
#' @param min_level minimum taxonomic level at which taxonomy
#' should be replaced. Default is the finest level (\code{"species"}).
#' @param force_species_name if \code{TRUE}, species not found in NCBI Taxonomy
#' will keep their original names instead of NAs.
#' @param verbose print information in the console.
#'
#' @return The reference database with the NCBI taxonomy
#' for the genus level and higher ranks.
#' (the original taxonomy above the genus level is removed).
#'
#' @examples
#' \dontrun{
#' lib <- read.csv(system.file("extdata", "baetidae_bold.csv", package = "refdb"))
#' lib <- refdb_set_fields_BOLD(lib)
#' refdb_set_ncbitax(lib)
#' }
#' @export
#'
#'
refdb_set_ncbitax <- function(x, min_level = "species", force_species_name = TRUE, verbose = TRUE) {

  check_fields(x, "taxonomy")
  ncbi_taxo <- ncbi_taxo_rank()
  valid_taxo <- valid_taxo_rank()
  x_taxo <- attributes(x)$refdb_fields$taxonomy
  x_taxo_sub_min_level <- x_taxo[!names(x_taxo) %in% valid_taxo[1:which(valid_taxo == min_level)]]
  x_taxo <- x_taxo[names(x_taxo) %in% valid_taxo[1:which(valid_taxo == min_level)]]
  x_taxo_sel <- x_taxo[stats::na.exclude(match(ncbi_taxo, names(x_taxo)))]
  x_taxo_sel <- rev(x_taxo_sel)

  # Columns matching a taxonomic field are suffixed
  colnames(x)[colnames(x) %in% x_taxo_sel] <-
    paste0(colnames(x)[colnames(x) %in% x_taxo_sel], "_original")

  x_taxo_sel_ori <- paste0(x_taxo_sel, "_original")
  names(x_taxo_sel_ori) <- names(x_taxo_sel)

  bt <- x[, x_taxo_sel_ori]
  bt <- dplyr::distinct(bt)

  # Collect NCBI IDs
  ncbi_ids <- vector("character", nrow(bt))
  for(i in seq_along(ncbi_ids)) {

    bt_i <- unlist(bt[i, ])
    ncbi_ids[i] <- NA
    lvl = 1

    while (is.na(ncbi_ids[i])) {

      if(is.na(bt_i[lvl])) {
        ncbi_ids[i] <- NA
      } else {
        if(verbose) cat("\rProcessing:", bt_i[lvl], rep(" ", 40))
        suppressMessages(
          uid <- taxize::get_uid_(bt_i[lvl], messages = FALSE)[[1]]$uid
        )

        # Deal with no match
        if(length(uid) == 0) {
          ncbi_ids[i] <- NA
        }

        # Deal with multiple matches
        if(length(uid > 1)) {
          tax <- get_ncbi_taxonomy(uid)
          score <- apply(tax[names(x_taxo_sel)], 1,
                         function(y) sum(y == bt_i, na.rm = TRUE)
          )
          ncbi_ids[i] <- tax$id[which.max(score)]
        }
      }
      lvl = lvl + 1
      # Deal with case where none of the taxonomic levels are found in NCBI
      # -> Use "unidentified"
      if(lvl > length(bt_i)) {
        ncbi_ids[i] <- 32644
        break
      }
    }
  }
  if(verbose) cat("\n\n")

  # Collect taxonomy for IDs
  bt_taxo <- get_ncbi_taxonomy(ncbi_ids)
  bt_taxo$id <- NULL
  taxo_field <- stats::setNames(colnames(bt_taxo), colnames(bt_taxo))

  # Columns of x which match a taxonomic name are removed
  x <- x[, !colnames(x) %in% colnames(bt_taxo)]

  # Expand results with left join
  bt_taxo <- dplyr::bind_cols(bt, bt_taxo)
  out <- dplyr::left_join(x, bt_taxo, by = unname(x_taxo_sel_ori))


 # Replace new NA values for species if original ones are not NA
  if(force_species_name) {
    sel <-
      !is.na(out[, x_taxo_sel_ori["species"], drop = TRUE]) &
      is.na(out[, taxo_field["species"], drop = TRUE])
    if(any(sel == TRUE)) {
      out[, taxo_field["species"]][sel, 1] <- out[, x_taxo_sel_ori["species"], drop = TRUE][sel]
    }
  }

  # Original taxonomy columns are removed
  out <- out[, !colnames(out) %in%
               c(x_taxo, x_taxo_sel_ori)]

  # Restore attributes and reorder taxonomic columns
  taxo_field <- c(taxo_field, x_taxo_sub_min_level)
  attributes(out)$refdb_fields$taxonomy <- taxo_field
  out <- out[, c(setdiff(colnames(out), taxo_field), taxo_field)]

  return(out)
}



#' Taxonomic ranks of the NCBI Taxonomy database
#'
#' @return a vector of ordered ranks
#'
ncbi_taxo_rank <- function() {
  c("superkingdom",
    "kingdom",
    "phylum",
    "subphylum",
    "class",
    "subclass",
    "infraclass",
    "order",
    "suborder",
    "infraorder",
    "superfamily",
    "family",
    "genus",
    "species")
}


#' Ranks considered as valid by refdb
#'
#' @return a vector of ordered ranks.
#' @export
#'
#' @examples
#' valid_taxo_rank()
#'
#' @references This is a simplified version of the
#' list \code{rank_ref} available in \pkg{taxize}.
#'
valid_taxo_rank <- function() {
  c(
    "domain",
    "superkingdom",
    "kingdom",
    "subkingdom",
    "infrakingdom",
    "phylum",
    "subphylum",
    "infradivision",
    "superclass",
    "class",
    "subclass",
    "infraclass",
    "subterclass",
    "parvclass",
    "megacohort",
    "supercohort",
    "cohort",
    "subcohort",
    "infracohort",
    "superorder",
    "order",
    "suborder",
    "infraorder",
    "parvorder",
    "superfamily",
    "family",
    "subfamily",
    "supertribe",
    "tribe",
    "subtribe",
    "genus",
    "subgenus",
    "section",
    "subsection",
    "series",
    "species",
    "infraspecies",
    "subspecies",
    "variety",
    "subvariety",
    "stirp",
    "form",
    "subform"
  )
}

