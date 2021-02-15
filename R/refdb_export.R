
#' Export reference database for DADA2
#'
#' Write reference database in formats which can be used with the
#' functions of the package \pkg{dada2}.
#'
#' @param x a reference database.
#' @param file a path to the file to be written.
#' @param mode character string to determine the type of file to produce.
#' Use \code{"taxonomy"} to produce a file for function \code{assignTaxonomy}
#' or \code{"species"} to produce a file for function \code{assignSpecies}.
#'
#' @export
#'
refdb_export_dada2 <- function(x, file, mode = "taxonomy") {
  check_fields(x, what = c("sequence", "taxonomy"))

  if(mode == "taxonomy") {
    col_tax <- attributes(x)$refdb_fields$taxonomy
    labs <- apply(x[, col_tax], 1, paste, collapse = ";")
    labs <- paste0(labs, ";")
  }

  if(mode == "species") {
    col_tax <- attributes(x)$refdb_fields$taxonomy["species"]
    col_id <- attributes(x)$refdb_fields$id

    labs <- apply(x[, c(col_id, col_tax)], 1, paste, collapse = " ")
    labs <- stringr::str_replace_all(labs, "(?<=[^ ]{1,1000} [^ ]{1,1000} [^ ]{1,1000}) ", "_")
    # Worst regex ever
  }

  seqs <- x[[attributes(x)$refdb_fields$sequence]]
  names(seqs) <- labs

  bioseq::write_fasta(seqs, file = file)

}
