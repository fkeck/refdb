
#' Compile a report with different checks
#'
#' This function produce an HTML report to investigate potential issues
#' in a reference database.
#'
#' @param x a reference database.
#' @param file the file (path) to write the report. If \code{NULL}
#' the report is written in the user temp directory.
#' @param view A logical. If \code{TRUE} (default), the file is instantly
#' opened in the web browser.
#'
#' @return
#' The function invisibly returns the file where the report was written.
#'
#' @export
#'
refdb_report <- function(x, file = NULL, view = TRUE) {

  if (!requireNamespace("DT", quietly = TRUE)) {
    stop("Package \"DT\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if(is.null(file)) {
    file <- tempfile("refdb_report_", fileext = ".html")
  }

  environment <- new.env()
  environment$x <- x

  template_path <- file.path(path.package("refdb"), "rmd_templates", "check_report.Rmd")

  rmarkdown::render(template_path,
                    output_format = "html_document",
                    output_file = file,
                    envir = environment)

  if(view) {
    utils::browseURL(file)
  }

  invisible(file)
}



#' Check for typos in taxonomic names
#'
#'This function uses the generalized Levenshtein (edit) distance
#'to identify possible issue with taxonomic names.
#'
#' @param x a reference database.
#' @param tol the edit distance below which two taxonomic names are reported.
#'
#' @return A list of two-columns tibbles reporting for each taxonomic level
#' the pairs of taxonomic names sharing the same upstream taxonomy and for
#' which the generalized Levenshtein (edit) distance is below
#' the \code{tol} value.
#' @export
#'
refdb_check_tax_typo <- function(x, tol = 1) {
  check_fields(x, "taxonomy")

  tax <- x[, attributes(x)$refdb$taxonomy]

  ff <- function(x, tol) {
    x <- unique(x)
    d <- utils::adist(x, x)
    d[upper.tri(d)] <- 0
    sel <- which(d > 0 & d <= tol, arr.ind = TRUE)
    sel <- apply(sel, 1, function(k) x[k])
    list(sel)
  }

  res <- vector("list", ncol(tax))

  for(i in seq_len(ncol(tax))) {
    tax <- dplyr::group_by(tax, !!!rlang::syms(colnames(tax)[seq(1, i-1)]))
    tax_ff <- dplyr::summarise(tax, tax_ff = ff(!!rlang::sym(colnames(tax)[i]), tol = tol))$tax_ff
    tax_ff <- tax_ff[sapply(tax_ff, length) > 0]
    if(length(tax_ff) == 0) {
      #res[[i]] <- NULL
    } else {
      tax_ff <- lapply(tax_ff, t)
      tax_ff <- do.call("rbind", tax_ff)
      colnames(tax_ff) <- c("Taxon 1", "Taxon 2")
      res[[i]] <- tibble::as_tibble(tax_ff)
    }
  }

  names(res) <- colnames(tax)

  return(res)

}





#' Check for conflicts in taxonomy
#'
#' @param x a reference database.
#'
#' @return A list of two-columns tibbles reporting for each taxonomic level
#' the taxa with identical names but different upstream taxonomy.
#'
#' @export
#'
refdb_check_tax_conflict <- function(x) {

  check_fields(x, "taxonomy")

  tax <- x[, attributes(x)$refdb$taxonomy]

  ff <- function(x) {
    x <- unique(x)
    list(x)
  }

  res <- vector("list", ncol(tax))

  for(i in seq_len(ncol(tax))) {
    tax <- dplyr::group_by(tax, !!rlang::sym(colnames(tax)[i]))
    tax_mod <- dplyr::mutate(tax,
                             refdb_pasted_cols = paste(!!!rlang::syms(colnames(tax)[seq(1, i-1)]),
                                                       sep = " > "))
    tax_ff <- dplyr::summarise(tax_mod, taxonomy = ff(.data$refdb_pasted_cols))
    tax_ff <- dplyr::filter(tax_ff, sapply(.data$taxonomy, length) > 1)
    colnames(tax_ff) <- c("Taxon", "Taxonomy")
    tax_ff <- dplyr::filter(tax_ff, !is.na(.data$Taxon))
    res[[i]] <- tidyr::unnest(tax_ff, .data$Taxonomy)
  }

  names(res) <- colnames(tax)

  return(res)

}



#' Check for conflicts in sequences
#'
#' @param x a reference database.
#' @param na_omit if \code{FALSE} conflicts involving
#' NA taxonomic names are also reported.
#'
#' @return A list of two-columns tibbles reporting duplicated
#' sequences with different taxonomy.
#'
#' @export
#'
refdb_check_seq_conflict <- function(x, na_omit = TRUE) {

  check_fields(x, what = c("taxonomy", "sequence"))

  tax_cols <- attributes(x)$refdb$taxonomy
  seq_col <- attributes(x)$refdb$sequence

  dat <- x[, c(tax_cols, seq_col)]

  dat <- dplyr::add_count(dat, !!rlang::sym(seq_col), name = "refdb_n_col")
  dat <- dplyr::filter(dat, .data$refdb_n_col > 1)
  dat$refdb_n_col <- NULL
  dat <- dplyr::group_by(dat, !!rlang::sym(seq_col))

  if(na_omit) {
    filter_tax <- function(x) {length(unique(x[!is.na(x)]))}
  } else {
    filter_tax <- function(x) {length(unique(x))}
  }

  dat <- tidyr::nest(dat)
  names(dat)[names(dat) == "data"] <- "dd"
  dat <- dplyr::mutate(dat,
                min_conf = sapply(.data$dd, function(x) {
                  s <- apply(x, 2, filter_tax)
                  names(tax_cols)[which(s > 1)][1]
                  }))
  dat <- dplyr::filter(dat, !is.na(.data$min_conf))
  dat <- tidyr::unnest(dat, cols = .data$dd)

  if(nrow(dat) == 0) {
    colnames(dat) <- c("sequence", "level_conflict", "taxonomy")
    return(dat)
  }

  dat <- dplyr::mutate(dat, refdb_pasted_cols = paste(!!!rlang::syms(tax_cols), sep = " > "))


  ff <- function(x) {
    list(unique(x))
  }

  dat <- dplyr::group_by(dat, !!rlang::sym(seq_col), .data$min_conf)
  dat <- dplyr::summarise(dat, taxonomy = ff(.data$refdb_pasted_cols))
  dat <- dplyr::filter(dat, sapply(.data$taxonomy, length) > 1)
  res <- tidyr::unnest(dat, .data$taxonomy)
  colnames(res) <- c("sequence", "level_conflict", "taxonomy")

  return(res)
}

