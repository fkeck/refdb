

#' Scores for filtering operations
#'
#' @param x a reference database
#'
#' @return a numeric vector
#'
.filter_seq_length <- function(x) {
  check_fields(x, "sequence")
  col <- attributes(x)$refdb$sequence
  bioseq::seq_nchar(x[[col]])
}

.filter_seq_ambiguous <- function(x, char) {
  check_fields(x, "sequence")
  col <- attributes(x)$refdb$sequence
  bioseq::seq_count_pattern(x[[col]], pattern = list(char))
}

.filter_seq_homopolymers <- function(x) {
  check_fields(x, "sequence")
  col <- attributes(x)$refdb$sequence
  res <- bioseq::seq_extract_pattern(x[[col]], pattern = "(.)\\1+")
  lapply(res, function(x) sort(unique(nchar(x))))
}

.filter_seq_duplicates <- function(x) {
  check_fields(x, what = c("sequence", "taxonomy"))
  col <- c(attributes(x)$refdb$sequence,
           attributes(x)$refdb$taxonomy)
  !duplicated.data.frame(x[, col])
}

.filter_seq_stopcodon <- function(x, code, codon_frame = NA) {
  check_fields(x, "sequence")
  col <- attributes(x)$refdb$sequence

  dna_seq <- x[[col]]
  dna_seq <- bioseq::seq_remove_pattern(dna_seq, pattern = "-")

  if(is.character(code)) {
    # If code is a column name
    code <- x[[code]]
  } else {
    #Otherwise it is assumed to be a numeric vector
    #indicating genetic code
    code <- rep_len(code, nrow(x))
  }

  if(!is.na(codon_frame)) {

    res <- mapply(function(x, y) {
      prot <- bioseq::seq_translate(x, code = y, codon_frame = codon_frame)
      bioseq::seq_count_pattern(prot, "\\*")
    },
    x = x[[col]],
    y = code,
    SIMPLIFY = TRUE,
    USE.NAMES = FALSE)

  } else {
    # We test the three frames
    stop_mat <- cbind(
      mapply(function(x, y) {
        prot <- bioseq::seq_translate(x, code = y, codon_frame = 1)
        bioseq::seq_count_pattern(prot, "\\*")
      },
      x = x[[col]],
      y = code,
      SIMPLIFY = TRUE,
      USE.NAMES = FALSE),
      mapply(function(x, y) {
        prot <- bioseq::seq_translate(x, code = y, codon_frame = 2)
        bioseq::seq_count_pattern(prot, "\\*")
      },
      x = x[[col]],
      y = code,
      SIMPLIFY = TRUE,
      USE.NAMES = FALSE),
      mapply(function(x, y) {
        prot <- bioseq::seq_translate(x, code = y, codon_frame = 3)
        bioseq::seq_count_pattern(prot, "\\*")
      },
      x = x[[col]],
      y = code,
      SIMPLIFY = TRUE,
      USE.NAMES = FALSE)
    )
    res <- apply(stop_mat, 1, min)
  }

  return(res)

}


.filter_tax_precision <- function(x) {
  check_fields(x, "taxonomy")
  col <- attributes(x)$refdb$taxonomy
  x_na <- refdb_clean_tax_NA(x)
  x_na <- x_na[, col]

  x_val <- lapply(seq_len(ncol(x_na)), function(x) {
    ifelse(is.na(x_na[, x]), 0, x)
  })

  apply(as.data.frame(x_val), 1, max)
}


.filter_ref_focus <- function(x) {
  check_fields(x, c("taxonomy", "reference"))
  col_tax <- attributes(x)$refdb$taxonomy
  col_ref <- attributes(x)$refdb$reference

  dat <- x[, c(col_ref, col_tax)]
  dat <- dplyr::group_by(dat, !!!rlang::syms(unname(col_ref)))

  dat <- dplyr::mutate(dat,
                       dplyr::across(.cols = dplyr::all_of(unname(col_tax)),
                                     .fns = dplyr::n_distinct))

  apply(dat[, col_tax], 1, function(x) sum(x > 1))
}

###############################


refdb_filter_seq_length <- function(x, min_len, max_len) {
  flt <- .filter_seq_length(x)
  sel <- flt >= min_len & flt <= max_len
  x[sel, ]
}



refdb_filter_seq_ambiguous <- function(x, max_ambig = 3, char = "N") {
  flt <- .filter_seq_ambiguous(x, char = char)
  sel <- flt <= max_ambig
  x[sel, ]
}



refdb_filter_seq_homopolymers <- function(x, max_len = 16) {
  flt <- .filter_seq_homopolymers(x)
  sel <- sapply(flt, max) <= max_len
  x[sel, ]
}


refdb_filter_seq_duplicates <- function(x) {
  sel <- .filter_seq_duplicates(x)
  x[sel, ]
}


# Can be renamed sample_
# Method random, max_genetic_div, central_genetic_div

# .filter_tax_subsample <- function(x, min_tax = "order", method, keep) {
#   check_fields(x, what = "taxonomy")
#   col <- attributes(x)$refdb$taxonomy
#   col <- col[seq(1, which(col == min_tax))]
#
#   if (method == "random") {
#     x[, col]
#
#   }
# }
#
# refdb_filter_seq_subsample <- function(x, method, keep) {
#   sel <- .filter_seq_subsample(x)
#   x[sel, ]
# }




refdb_filter_seq_stopcodon <- function(x, max_stop = 0, code, codon_frame = NA) {
  flt <- .filter_seq_stopcodon(x, code = code, codon_frame = codon_frame)
  sel <- flt <= max_stop
  x[sel, ]
}


refdb_filter_tax_precision <- function(x, min_tax) {
  flt <- .filter_tax_precision(x)

  min_tax <- which(min_tax == attributes(x)$refdb$taxonomy)
  sel <- flt >= min_tax

  x[sel, ]
}


refdb_filter_ref_focus <- function(x, min_tax) {
  flt <- .filter_ref_precision(x)

  min_tax <- which(min_tax == attributes(x)$refdb$taxonomy)
  sel <- flt >= min_tax

  x[sel, ]
}

# refdb_filter_seq_primer()

# #Require alignment:
# refdb_filter_seq_dist(x, max_dist)
#
# #Require alignment or phylogenetic tree
# refdb_filter_seq_phylo(x, phylo, max_dist)
#
# #Require external RDP software
# refdb_filter_seq_selfassign(x, max_dist, exec)
