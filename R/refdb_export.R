
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

  if(mode == "taxonomy") {
    check_fields(x, what = c("sequence", "taxonomy"))
    col_tax <- attributes(x)$refdb_fields$taxonomy
    labs <- apply(x[, col_tax], 1, paste, collapse = ";")
    labs <- paste0(labs, ";")
  }

  if(mode == "species") {
    check_fields(x, what = c("sequence", "taxonomy", "id"))
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





#' Export reference database for Mothur
#'
#' Write a reference database in formats which can be used
#' with \code{Mothur}..
#'
#' @param x a reference database.
#' @param file a file path. This will be used to create
#'  a .fasta file and a .txt file.
#'
#' @export
#'
refdb_export_mothur <- function(x, file) {
  check_fields(x, what = c("sequence", "taxonomy", "id"))

  col_tax <- attributes(x)$refdb_fields$taxonomy
  col_id <- attributes(x)$refdb_fields$id
  col_seq <- attributes(x)$refdb_fields$sequence

  labs <- x[[col_id]]
  labs <- stringr::str_replace_all(labs, "[:blank:]", "_")

  seqs <- x[[col_seq]]
  names(seqs) <- labs

  tax <- apply(x[, col_tax], 1, function(x) stringr::str_replace_all(paste(x, collapse = ";"), "[:blank:]", "_"))
  tax <- paste0(labs, "\t", tax, ";")


  file_fas <- paste0(file, ".fasta")
  file_txt <- paste0(file, ".txt")

  bioseq::write_fasta(seqs, file = file_fas)
  readr::write_lines(tax, file = file_txt)
}


