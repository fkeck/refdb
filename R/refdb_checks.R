
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
#' @examples
#'
#' \dontrun{
#' lib <- read.csv(system.file("extdata", "ephem.csv", package = "refdb"))
#' lib <- refdb_set_fields(lib,
#'                         taxonomy = c(family = "family_name",
#'                             genus = "genus_name",
#'                             species = "species_name"),
#'                         sequence = "DNA_seq",
#'                         marker = "marker")
#' tmp <- tempfile()
#'
#' refdb_report(lib, tmp, view = FALSE)
#' }
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
#'
#' @examples
#' lib <- read.csv(system.file("extdata", "ephem.csv", package = "refdb"))
#' lib <- refdb_set_fields(lib,
#'                         taxonomy = c(family = "family_name",
#'                             genus = "genus_name",
#'                             species = "species_name"),
#'                         sequence = "DNA_seq",
#'                         marker = "marker")
#' refdb_check_tax_typo(lib)
#'
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
#' @examples
#' lib <- read.csv(system.file("extdata", "ephem.csv", package = "refdb"))
#' lib <- refdb_set_fields(lib,
#'                         taxonomy = c(family = "family_name",
#'                             genus = "genus_name",
#'                             species = "species_name"),
#'                         sequence = "DNA_seq",
#'                         marker = "marker")
#' refdb_check_tax_conflict(lib)
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
#' @examples
#' lib <- read.csv(system.file("extdata", "ephem.csv", package = "refdb"))
#' lib <- refdb_set_fields(lib,
#'                         taxonomy = c(family = "family_name",
#'                             genus = "genus_name",
#'                             species = "species_name"),
#'                         sequence = "DNA_seq",
#'                         marker = "marker")
#' refdb_check_seq_conflict(lib)
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





#' Check for genetic homogeneity of taxa
#'
#' This function assesses the genetic similarity among
#' sequences within each taxa.
#' It takes user defined thresholds (one threshold per taxonomic level)
#' to warn about sequences which are singularly different
#' (based on median distance) from the others.
#' Sequences in the reference database must be aligned.
#'
#' @param x a reference database (sequences must be aligned).
#' @param levels a named vector of genetic similarity thresholds.
#' Names must correspond to taxonomic levels (taxonomic fields)
#' and values must be included in the interval [0, 1].
#' For example to assess homogeneity at 5 percents (within species) and
#' 10 percents (within genus): \code{levels = c(species = 0.05, genus = 0.1)}
#' @param min_n_seq the minimum number of sequences for a
#' taxon to be tested.
#'
#' @details
#' For every tested taxonomic levels, the algorithm
#' checks all sequences in every taxa
#' (for which the total number of sequence is > \code{min_n_seq})
#' In each taxon, the pairwise distance matrix among all the sequences
#' belonging to this taxon is computed. A sequence is tagged as suspicious
#' and returned by the function
#' if its median genetic distance from the other sequences is higher than
#' the threshold set by the user (\code{levels} argument).
#'
#'
#' @return
#' A dataframe reporting suspicious sequences whose median distance
#' to other sequences of the same taxon is greater than
#' the specified threshold. The first column "level_threshold_homogeneity"
#' indicates the lowest taxonomic level for which the threshold has been
#' exceeded and the second column "value_threshold_homogeneity" gives
#' the computed median distance.
#'
#' @examples
#' lib <- read.csv(system.file("extdata", "homogeneity.csv", package = "refdb"))
#' lib <- refdb_set_fields_BOLD(lib)
#' refdb_check_seq_homogeneity(lib, levels = c(species = 0.05, genus = 0.1))
#'
#' @export
#'
refdb_check_seq_homogeneity <- function(x, levels, min_n_seq = 3) {

  check_fields(x, what = c("taxonomy", "sequence"))

  levels <- levels[order(match(names(levels), rev(ncbi_taxo_rank())))]

  tax_cols <- attributes(x)$refdb$taxonomy
  seq_col <- attributes(x)$refdb$sequence

  # A weak test to know if sequences are aligned
  if(!all(nchar(x[, seq_col]) == nchar(x[1, seq_col]))) {
    stop("Sequences must be aligned.")
  }

  dat <- x[, c(tax_cols, seq_col)]

  median_dist <- function(x, min_n_seq) {
    if(length(x) < min_n_seq) {
      return(NA)
    }
    seq_i_g <- bioseq::as_DNAbin(x)
    seq_dist_i_g <- ape::dist.dna(seq_i_g, model = "raw", as.matrix = TRUE)
    unname(apply(seq_dist_i_g, 2, stats::median))
  }

  res <- vector("list", length(levels))
  names(res) <- names(levels)

  for (i in names(levels)){
    tax_i <- 1:which(names(tax_cols) == i)
    dat_i <- dat[, c(tax_i, ncol(dat))]
    dat_i <- dplyr::group_by(dat_i, dplyr::across(1:max(tax_i)))
    dat_i <- dplyr::mutate(dat_i, med_dist = median_dist(!!rlang::sym(seq_col), min_n_seq = min_n_seq))
    res[[i]] <- dat_i$med_dist
  }

  sel <- mapply(function(x, y) x > y, x = res, y = levels)
  level_threshold_homogeneity <- apply(sel, 1, function(x) colnames(sel)[x][1])
  value_threshold_homogeneity <- apply(as.data.frame(res), 1, function(x) x[x > levels][1])
  res <- cbind(level_threshold_homogeneity, value_threshold_homogeneity, x)
  res <- dplyr::filter(res, !is.na(level_threshold_homogeneity))
  return(res)
}
