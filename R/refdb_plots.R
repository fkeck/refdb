

#' Create a graph from a taxonomic table
#'
#'Create a graph representation from a taxonomic classification
#'included in a reference
#'database. For this function to work, taxonomic fields must be set.
#'
#'
#' @param x a reference database (tibble).
#' @param cols an optional vector of column names to use a subset of columns.
#'
#' @return
#' An \pkg{igraph} object representing taxonomic relationships.
#'
#'
igraph_from_taxo <- function(x, cols = NULL) {

  check_fields(x, what = "taxonomy")
  col_tax <-attributes(x)$refdb_fields$taxonomy
  x <- tibble::as_tibble(x)

  x <- x[, col_tax]

  if(!is.null(cols)) {
    cols <- cols[order(match(cols, col_tax))]
    x <- x[, cols]
  }

  if(length(unique(unlist(x[, 1]))) > 1) {
    x <- tibble::tibble(Root = "Root", x)
  }

  res <- x
  for(i in 1:ncol(res)) {
    res[, i] <- do.call("paste", c(x[, 1:i], sep = ">"))
  }

  # Nodes
  res_tidy <- vector("list", ncol(res))
  for(i in 1:ncol(res)) {
    res_tidy[[i]] <- dplyr::group_by(res[, i], dplyr::across())
    res_tidy[[i]] <- dplyr::count(res_tidy[[i]])
    res_tidy[[i]] <- dplyr::ungroup(res_tidy[[i]])
    names(res_tidy[[i]]) <- c("taxonomy_path", "n")
    res_tidy[[i]]$taxonomy_rank <- names(res)[i]
  }
  res_tidy <- do.call("rbind", res_tidy)
  res_tidy$taxonomy_depth <- stringr::str_count(res_tidy$taxonomy_path, ">")
  res_tidy$terminal <- stringr::str_remove(res_tidy$taxonomy_path, ".+>")
  res_tidy <- dplyr::group_by(res_tidy, .data$taxonomy_rank)
  res_tidy <- dplyr::mutate(res_tidy, freq_by_rank = .data$n/sum(.data$n))
  res_tidy <- dplyr::ungroup(res_tidy)


  # Edges
  res_derep <- dplyr::distinct(res)
  graph_def_edges <- apply(res_derep, 1, function(r) {
    r <- r[!stringr::str_detect(r, ">NA$")]
    if(length(r) < 2) return(NULL)
    edges <- matrix(c(r[1:(length(r)-1)], r[2:(length(r))]),
                    ncol = 2, dimnames = list(NULL, c("from", "to")))
    as.data.frame(edges)
  })
  graph_def_edges <- do.call(rbind, graph_def_edges)
  graph_def_edges <- dplyr::distinct(graph_def_edges)
  graph_def_edges$branch_length <-
    stringr::str_count(graph_def_edges$to, ">") -
    stringr::str_count(graph_def_edges$from, ">")

  res_tidy <- res_tidy[res_tidy$taxonomy_path %in%
                         c(graph_def_edges$from, graph_def_edges$to), ]

  g <- igraph::graph_from_data_frame(graph_def_edges,
                                     vertices = res_tidy)

  return(g)
}





#' Reference database treemap
#'
#' Represent the hierarchical structure of the taxonomic information
#' of a reference database as a set of nested rectangles (treemap).
#'
#' @param x a reference database.
#' @param cols a vector of column names referring to taxonomic levels
#' to include in the treemap. If not provided (\code{NULL})
#' the function tries to find a relevant subset of columns.
#' @param freq_labels a numeric vector of length two to adjust
#' the number of printed labels (see Details).
#'
#' @details
#' TODO
#'
#' Only the columns provided in the \code{cols} argument are represented
#' in the treemap. Large labels are printed for the highest rank, while light
#' text labels are printed for the lowest rank
#' Intermediate ranks are drawn but their names are not shown.
#'
#' The number of labels printed are determined by \code{freq_labels}.
#' The first value gives the threshold
#' for the highest rank (large labels) and the second value gives
#' the threshold for the lowest rank (light text labels).
#'
#' The underlying graph is computed using the non-exported function
#' \code{igraph_from_taxo}.
#'
#' @return
#' A \pkg{ggplot2} (\pkg{ggraph}) object. This means the plot
#' can be further customized using \pkg{ggplot2} compatible functions.
#' @export
#'
refdb_plot_tax_treemap <- function(x, cols = NULL, freq_labels = c(0.01, 0.003)) {

  check_fields(x, what = "taxonomy")
  col_tax <- attributes(x)$refdb_fields$taxonomy


  if(is.null(cols)) {
    cols_n <- apply(x[, col_tax], 2, function(x) length(unique(x)))
    cols <- names(cols_n[cols_n > 50 & cols_n < 5000])
    if(length(cols) < 2) {
      cols <- names(cols_n)[c(which.max(cols_n) - 1, which.max(cols_n))]
      }
    cat("Selected columns: ", cols, "\n")
  }

  g <- igraph_from_taxo(x, cols = cols)

# Get labels for the nth taxonomic level
# from taxonomic string separated by ">"
# Root is zero
  taxo_nth_level <- function(x, level){
    out <- stringr::str_split_fixed(x, ">", n = Inf)
    out[out == "" | out == "NA"] <- NA
    out[, level + 1]
  }

  ggraph::ggraph(g, 'treemap', weight =  .data$n) +
    ggraph::geom_node_tile(ggplot2::aes(fill = taxo_nth_level(.data$name, 1),
                                        alpha = .data$n),
                           size = 0, show.legend = FALSE) +
    ggraph::geom_node_tile(ggplot2::aes(size = .data$taxonomy_depth),
                           colour = 'white') +
    ggraph::geom_node_text(ggplot2::aes(label = ifelse(.data$taxonomy_depth == max(.data$taxonomy_depth) &
                                                         .data$freq_by_rank > freq_labels[2],
                                                       .data$terminal,
                                                       NA_character_)),
                           size = 3, color = "white", alpha = 0.8) +
    ggraph::geom_node_label(ggplot2::aes(label = ifelse(.data$taxonomy_depth == 1 &
                                                          .data$freq_by_rank > freq_labels[1],
                                                        .data$terminal,
                                                        NA_character_))) +
    ggplot2::scale_alpha(range = c(1, 0.5), guide = 'none') +
    ggplot2::scale_size(range = c(4, 0.2), guide = 'none') +
    ggplot2::theme_void()

}





#' Reference database taxonomy tree
#'
#' Represent the hierarchical structure of the taxonomic information
#' of a reference database as a tree.
#'
#' @param x a reference database.
#' @param leaf_col a column name referring to the taxonomic level
#' for the leaves of the tree. If not provided (\code{NULL})
#' the function tries to find a relevant level.
#' @param color_col a column name referring to the taxonomic level
#' for the color of the leaves (must be higher or equal to the level
#' of \code{leaf_col}). If not provided (\code{NULL})
#' the function tries to find a relevant level.
#' @param freq_labels a numeric value to adjust
#' the number of printed labels (minimum frequency).
#' Default is zero which means all non-NA labels are printed.
#' @param expand_plot a value to expand the limits of the plot.
#' Useful if the labels are too long.
#'
#' @details
#' The underlying graph is computed using the non-exported function
#' \code{igraph_from_taxo}.
#'
#' @return
#' A \pkg{ggplot2} (\pkg{ggraph}) object. This means the plot
#' can be further customized using \pkg{ggplot2} compatible functions.
#' @export
#'
refdb_plot_tax_tree <- function(x,
                                leaf_col = NULL,
                                color_col = NULL,
                                freq_labels = 0,
                                expand_plot = 0.5) {

  check_fields(x, what = "taxonomy")
  col_tax <- attributes(x)$refdb_fields$taxonomy

  # Get labels for the nth taxonomic level
  # from taxonomic string separated by ">"
  # Root is zero
  taxo_nth_level <- function(x, level){
    out <- stringr::str_split_fixed(x, ">", n = Inf)
    out[out == "" | out == "NA"] <- NA
    if(all(out[, 1] == "Root")) {
      out <- out[, level + 1]
    } else {
      out <- out[, level]
    }
    return(out)
  }


  if(is.null(leaf_col)) {
    leaf_col <- apply(x[, col_tax], 2, function(x) length(unique(x)))
    leaf_col <- names(leaf_col[leaf_col < 200])
    cat("Selected rank columns for the tree: ", leaf_col, "\n")
  } else {
    leaf_col <- col_tax[1:which(col_tax == leaf_col)]
  }

  if(is.null(color_col)) {
    color_col <- apply(x[, leaf_col], 2, function(x) length(unique(x)))
    color_col <- names(color_col[color_col < 25])
    color_col <- color_col[length(color_col)]
    cat("Selected rank column for the color: ", color_col, "\n")
  } else {
    color_col <- col_tax[col_tax == color_col]
  }

  g <- igraph_from_taxo(x, cols = leaf_col)

  lvl_color_col <- ifelse(length(color_col) == 0, 0, which(col_tax == color_col))
  igraph::V(g)$leaf_group <- taxo_nth_level(igraph::V(g)$name, lvl_color_col)

  ggraph::ggraph(g, 'dendrogram', circular = TRUE) +
    ggraph::geom_edge_diagonal(colour = "grey30") +
    ggraph::geom_node_point(ggplot2::aes(filter = .data$leaf,
                                         size = .data$n,
                                         color = .data$leaf_group),
                            alpha = 0.5) +
    ggraph::geom_node_text(ggplot2::aes(x = .data$x*1.05,
                                        y = .data$y*1.05,
                                        label = ifelse(.data$taxonomy_depth == max(.data$taxonomy_depth) &
                                                         .data$freq_by_rank > freq_labels,
                                                       .data$terminal,
                                                       NA_character_),
                                        angle = -((- ggraph::node_angle(.data$x, .data$y) + 90) %% 180) + 90),
                           size = 3,
                           color = "black",
                           hjust='outward') +
    ggplot2::scale_size_continuous(range = c(1, 8) ) +
    ggplot2::coord_fixed() +
    ggplot2::expand_limits(x = c(-1 - expand_plot, 1 + expand_plot),
                           y = c(-1 - expand_plot, 1 + expand_plot)) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "none")

}




#' Plot an histogram of sequence lengths
#'
#' @param x a reference database
#' @param remove_gaps a logical (default \code{TRUE} to control whether
#' gaps (-) should be removed prior computing sequence lengths.
#'
#' @return
#' A ggplot object. This means the plot
#' can be further customized using \pkg{ggplot2} compatible functions.
#'
#' @export
#'
refdb_plot_seqlen_hist <- function(x, remove_gaps = TRUE) {
  check_fields(x, what = "sequence")
  col_seq <- attributes(x)$refdb_fields$sequence

  seqs <- x[[col_seq]]

  if(remove_gaps) {
    seqs <- bioseq::seq_remove_pattern(seqs, "-")
  }

  dat_len <- data.frame(Length = bioseq::seq_nchar(seqs))

  ggplot2::ggplot(dat_len) +
    ggplot2::geom_histogram(ggplot2::aes(.data$Length)) +
    ggplot2::xlab("Sequence length (bp)") +
    ggplot2::ylab("Count") +
    ggplot2::theme_bw()

}





