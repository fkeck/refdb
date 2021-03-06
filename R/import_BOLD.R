

#' Download and import BOLD records
#'
#' This function allows to search and download data
#' from the the BOLD database.
#' Additionally it allows to convert the BOLD taxonomic classification
#' into the NCBI classification using the NCBI Taxonomy database.
#'
#' @param taxon,ids,bin,container,institutions,researchers,geo terms
#' used to search the BOLD database. See \link[bold]{bold_seqspec} from
#' which they are inherited.
#' @param ncbi_taxo if \code{TRUE} (the default) the BOLD taxonomy is
#' replaced by the NCBI taxonomy using \code{\link{refdb_ncbi_taxonomy}}.
#' @param full if \code{TRUE} the complete data are returned for every records.
#' If \code{FALSE} (the default) only a subset of the most
#' important columns is returned.
#'
#' @details
#' This function uses several functions of the \pkg{bold}
#' package to interface with the BOLD API. The NCBI Taxonomy
#' database is queried using the \pkg{rentrez} package.
#'
#' @section Large requests:
#' As detailed in the manual of bold_seqspec, some large requests
#' (e.g. high ranked taxa) can lead to errors. In that case one
#' strategy can be to sequentially download data for lower rank
#' taxa. See \url{https://docs.ropensci.org/bold/#large-data}.
#'
#' @return A tibble.
#'
#' @seealso link[bold]{bold_stats} and \link[bold]{bold_seqspec}
#' from \pkg{bold} which are used internally.
#'
#' @export
#'
#'
#'
refdb_import_BOLD <- function(taxon = NULL,
                              ids = NULL,
                              bin = NULL,
                              container = NULL,
                              institutions = NULL,
                              researchers = NULL,
                              geo = NULL,
                              ncbi_taxo = TRUE,
                              full = FALSE) {

  query_stats <- bold::bold_stats(taxon = taxon,
                                  ids = ids,
                                  bin = bin,
                                  container = container,
                                  institutions = institutions,
                                  researchers = researchers,
                                  geo = geo)

  if(query_stats$total_records == 0) {
    cat("No sequence found\n")
    return(NULL)
  }

  cat("Downloading", query_stats$total_records, "sequences from BOLD...\n")

  tryCatch(
    recs <- bold::bold_seqspec(taxon = taxon,
                               ids = ids,
                               bin = bin,
                               container = container,
                               institutions = institutions,
                               researchers = researchers,
                               geo = geo),
    error = function(c) {
      stop(
      "\nAn error occured while trying to download data from BOLD servers\n",
      "For large requests check the manual (?refdb_import_BOLD).\n"
      )
    }
  )
  out <- tibble::as_tibble(recs)


  out <- tibble::tibble(source = "BOLD", out)

  if (full == FALSE) {
    out <- dplyr::select(out,
                         .data$source,
                         .data$sequenceID,
                         .data$markercode,
                         .data$phylum_name,
                         .data$class_name,
                         .data$order_name,
                         .data$family_name,
                         .data$subfamily_name,
                         .data$genus_name,
                         .data$species_name,
                         .data$subspecies_name, # Correspond to NCBI organism
                         .data$nucleotides,
                         .data$country,
                         .data$province_state,
                         .data$lat,
                         .data$lon)
  }

  # Empty strings as NA
  out[out == ""] <- NA

  out <- refdb_set_fields(out, source = "source",
                          id = "sequenceID",
                          organism = "species_name",
                          taxonomy = c(phylum = "phylum_name",
                                       class = "class_name",
                                       order = "order_name",
                                       family = "family_name",
                                       subfamily = "subfamily_name",
                                       genus = "genus_name"),
                          sequence = "nucleotides",
                          marker = "markercode")

  if (ncbi_taxo) {
    out <- refdb_ncbi_taxonomy(out)
  }

  return(out)
}


# Attempt to implement data collection recursively in the taxonomy
# Not working yet!
# bold_seqspec_large <- function(taxon, ncbi_taxo, full) {
#
#   query_stats <- bold::bold_stats(taxon = taxon)
#
#   tax_levels <- c("order", "family", "genus", "species")
#   top_tax <- sapply(query_stats[tax_levels], function(x) x$count)
#
#   if(top_tax["order"] > 1) {
#     down_to_tax <- "order"
#   } else {
#     top_tax <- dplyr::last(names(which(top_tax == min(top_tax))))
#     down_to_tax <- tax_levels[which(tax_levels == top_tax) + 1]
#   }
#
#   suppressMessages(
#     tax_ds <- taxize::downstream(taxon, db = "bold", downto = "genus")
#   )
#
#   out <- lapply(tax_ds[[taxon]]$name, refdb_import_BOLD,
#                 ncbi_taxo = ncbi_taxo,
#                 full = full)
#   dplyr::bind_rows(out)
# }
