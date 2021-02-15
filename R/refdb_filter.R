

#' Scores for filtering operations
#'
#' @param x a reference database
#'
#' @return a numeric vector
#' @name filter_scores
NULL


#' @rdname filter_scores
.filter_seq_length <- function(x) {
  check_fields(x, "sequence")
  col <- attributes(x)$refdb$sequence
  bioseq::seq_nchar(x[[col]])
}

.filter_seq_ambiguous <- function(x, char) {
  check_fields(x, "sequence")
  col <- attributes(x)$refdb$sequence
  bioseq::seq_count_pattern(x[[col]], pattern = list(char))
} # Returns number of ambiguous char

.filter_seq_homopolymers <- function(x) {
  check_fields(x, "sequence")
  col <- attributes(x)$refdb$sequence
  res <- bioseq::seq_extract_pattern(x[[col]], pattern = "(.)\\1+")
  sapply(res, function(x) max(sort(unique(nchar(x)))))
} # Returns length of longest homopolymer

.filter_seq_duplicates <- function(x) {
  check_fields(x, what = c("sequence", "taxonomy"))
  col <- c(attributes(x)$refdb$sequence,
           attributes(x)$refdb$taxonomy)
  !duplicated.data.frame(x[, col])
} # Return logical vector: FALSE are duplicates

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
} # Return number of stop codons



.filter_seq_primer <- function(x, primer) {
  check_fields(x, "sequence")
  col <- attributes(x)$refdb$sequence

  dna_seq <- x[[col]]
  dna_seq <- bioseq::seq_remove_pattern(dna_seq, pattern = "-")

  strd <- bioseq:::seq_afind(dna_seq, bioseq::as_dna(primer))

  res <- apply(strd$distance, 1, min) /
    nchar(apply(cbind(1:nrow(strd$distance),
                      apply(strd$distance, 1, which.min)),
                1, function(x) strd$match[x[[1]], x[[2]]]))
  return(res)
} # Return distance to primer (frequency of character [0-1])



.filter_tax_precision <- function(x) {
  check_fields(x, "taxonomy")
  col <- attributes(x)$refdb$taxonomy
  x_na <- refdb_clean_tax_NA(x)
  x_na <- x_na[, col]

  x_val <- lapply(seq_len(ncol(x_na)), function(x) {
    ifelse(is.na(x_na[, x]), 0, x)
  })

  apply(as.data.frame(x_val), 1, max)
} # Return taxonomic precision (number of last taxonomic column)
  # High value indicate more precise identification

.filter_ref_scope <- function(x) {
  check_fields(x, c("taxonomy", "reference"))
  col_tax <- attributes(x)$refdb$taxonomy
  col_ref <- attributes(x)$refdb$reference

  dat <- x[, c(col_ref, col_tax)]
  dat <- dplyr::group_by(dat, !!!rlang::syms(unname(col_ref)))

  dat <- dplyr::mutate(dat,
                       dplyr::across(.cols = dplyr::all_of(unname(col_tax)),
                                     .fns = dplyr::n_distinct))

  res <- apply(dat[, col_tax], 1, function(x) sum(x > 1))
  res[is.na(dat[[col_ref]])] <- NA
  length(col_tax) - res
} # Return the study taxonomic scope (number of last taxonomic column)
# High value indicate more narrow studies
# Returns NA if reference is NA



.filter_seq_dist <- function(x) {
  check_fields(x, c("taxonomy", "sequence"))
  x <- refdb_clean_tax_NA(x)
  x_tax <- x[, attributes(x)$refdb_fields$taxonomy]
  x_seq <- x[, attributes(x)$refdb_fields$sequence, drop = TRUE]
  x_tax_precision <- .filter_tax_precision(x)
  g <- igraph_from_taxo(x)

  seq_in_out <- cbind(stringr::str_locate(x_seq, "^[-N]+")[, 2],
                      stringr::str_locate(x_seq, "[-N]+$")[, 1])
  seq_in_out[, 1][is.na(seq_in_out[, 1])] <- 1
  seq_in_out[, 2][is.na(seq_in_out[, 2])] <- nchar(x_seq[1])

  out <- vector(mode = "logical", length = length(x_seq))

  for(i in seq_along(x_seq)) {
    i_seq <- x_seq[i]
    i_seq_in_out <- seq_in_out[i, ]
    i_tax <- unlist(x_tax[i, ])
    i_tax_precision <- x_tax_precision[i]

    if(length(unique(unlist(x_tax[, 1]))) > 1) {
      i_tax_path <- c("Root", i_tax)
    } else {
      i_tax_path <- i_tax
    }
    i_tax_path <- paste(i_tax_path, collapse = ">")
    i_tax_path <- stringr::str_remove(i_tax_path, pattern = "(>NA)+$")


    sample_idx <- intersect(
      which(seq_in_out[, 1] <= i_seq_in_out[1]),
      which(seq_in_out[, 2] >= i_seq_in_out[2])
    )
    sample_idx <- intersect(
      sample_idx,
      which(x_tax_precision == i_tax_precision)
    )
    if(length(sample_idx) < 100 | is.na(x_seq[i])) {
      out[i] <- 0
      next
    }
    sample_idx <- sample(sample_idx, size = 100)
    sample_tax_path <- x_tax[sample_idx, ]
    if(length(unique(unlist(x_tax[, 1]))) > 1) {
      sample_tax_path <- tibble::tibble(Root = "Root", sample_tax_path)
    }

    sample_tax_path <- apply(x_tax[sample_idx, ], 1,
                             function(x) {
                               res <- paste(x, collapse = ">")
                               res <- stringr::str_remove(res, pattern = "(>NA)+$")
                               })

    sample_seq <- x_seq[sample_idx]

    tax_dist_i <- igraph::distances(graph = g,
                                    v = which(igraph::V(g)$name == i_tax_path),
                                    to = which(igraph::V(g)$name %in%  sample_tax_path))

    tax_dist_i <- tax_dist_i[ , sample_tax_path, drop = TRUE]

    seq_dist_i <- ape::dist.dna(bioseq::as_DNAbin(c(i_seq, sample_seq)),
                                model = "raw", as.matrix = TRUE)[-1, 1]

    sel <- seq_dist_i < quantile(seq_dist_i, probs = 0.95)

    seq_dist_i <- seq_dist_i[sel]
    tax_dist_i <- tax_dist_i[sel]

    median_seq_dist <- tapply(seq_dist_i, tax_dist_i, median)
    out[i] <- sum(median_seq_dist[1] > median_seq_dist[-1])

    if(all(tax_dist_i == tax_dist_i[1]) | all(seq_dist_i == seq_dist_i[1])) {
      out[i] <- 0
      next
    }

    if(cor.test(tax_dist_i, seq_dist_i)$conf.int[1] > 0) {
      out[i] <- 0
    }

    if(out[i] > 2) plot(tax_dist_i, seq_dist_i, main = i)

    cat("\rProcessing sequences: ", i, " (", floor(i/length(x_seq)*100), "%)", rep(" ", 40), sep = "")
  }

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
  sel <- flt <= max_len
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



refdb_filter_seq_primer <- function(x, primer_forward = NULL,
                                    primer_reverse = NULL,
                                    max_error_forward = 0.1,
                                    max_error_reverse = 0.1) {
  sel <- rep(TRUE, nrow(x))

  if(!is.null(primer_forward)) {
    flt_fwd <- .filter_seq_primer(x, primer = primer_forward)
    sel <- sel & (flt_fwd <= max_error_forward)
  }

  if(!is.null(primer_reverse)) {
    flt_rev <- .filter_seq_primer(x, primer = primer_reverse)
    sel <- sel & (flt_rev <= max_error_reverse)
  }

  x[sel, ]
}



refdb_filter_tax_precision <- function(x, min_tax) {
  flt <- .filter_tax_precision(x)

  min_tax <- which(min_tax == attributes(x)$refdb$taxonomy)
  sel <- flt >= min_tax
  x[sel, ]
}


#' Filter records by taxonomic scope of studies
#'
#' @param x a reference database (tibble).
#' @param max_tax the maximum taxonomic focus of the study.
#'
#' @details
#' A reference field (one ore more columns) must be set to use
#' this function. If reference is not available (NA) for a record,
#' the record is not dropped.
#'
#' @return
#' a reference database (tibble).
#'
#' @export
#'
refdb_filter_ref_scope <- function(x, max_tax) {
  flt <- .filter_ref_scope(x)

  min_tax <- which(min_tax == attributes(x)$refdb$taxonomy)
  sel <- flt >= min_tax
  sel[is.na(sel)] <- TRUE
  x[sel, ]
}


refdb_filter_seq_dist <- function(x, max_dist){

}



#

# #Require alignment:
# refdb_filter_seq_dist(x, max_dist)
#
# #Require alignment or phylogenetic tree
# refdb_filter_seq_phylo(x, phylo, max_dist)
#
# #Require external RDP software
# refdb_filter_seq_selfassign(x, max_dist, exec)




