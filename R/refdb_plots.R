

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
refdb_plot_tax_treemap <- function(x, cols = NULL,
                                   freq_labels = c(0.01, 0.003)) {

  check_fields(x, what = "taxonomy")
  col_tax <- attributes(x)$refdb_fields$taxonomy


  if(is.null(cols)) {
    cols_n <- apply(x[, col_tax], 2, function(x) length(unique(x)))
    cols <- names(cols_n[cols_n > 50 & cols_n < 5000])
    if(length(cols) < 2) {
      cols <- names(cols_n)[c(which.max(cols_n) - 1, which.max(cols_n))]
    }
      message("Selected columns: ", paste(cols, collapse = ", "))
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
    message("Selected rank columns for the tree: ", paste(leaf_col, collapse = ", "))
  } else {
    leaf_col <- col_tax[1:which(col_tax == leaf_col)]
  }

  if(is.null(color_col)) {
    color_col <- apply(x[, leaf_col], 2, function(x) length(unique(x)))
    color_col <- names(color_col[color_col < 25])
    color_col <- color_col[length(color_col)]
    message("Selected rank column for the color: ", paste(color_col, collapse = ", "))
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




#' Plot an interactive map
#'
#' This functions generate an interactive maps showing
#' the location of the records of a reference database.
#' Note that only records with
#' latitude and longitude data will be displayed.
#'
#' @param x a reference database.
#'
#' @return
#' An interactive map object from the \pkg{leaflet} package.
#'
#' @export
#'
refdb_plot_map <- function(x) {

  check_fields(x, c("taxonomy", "sequence", "latitude", "longitude"))

  if(!"id" %in% names(attributes(x)$refdb)) {
    x$id <- seq_len(nrow(x))
    x <- refdb_set_fields(x, id = "id")
  }

  check_fields(x, "id")

  tax_cols <- attributes(x)$refdb$taxonomy
  id_col <- attributes(x)$refdb$id
  seq_col <- attributes(x)$refdb$sequence
  lat_col <- attributes(x)$refdb$latitude
  lon_col <- attributes(x)$refdb$longitude

  tax_last <- .filter_tax_precision(x)
  tax <- x[tax_cols]

  tax_last <- mapply(function(x, y) tax[y, x], x = tax_last, y = seq_len(nrow(tax)))
  tax_last <- unlist(tax_last)

  popup <- paste0("<h3>ID: ", x[[id_col]], "</h3>",
                  "<b>", stringr::str_to_title(names(tax_last)),
                  ": </b>", tax_last,
                  "<br><b>Sequence length: </b>", nchar(x[[seq_col]]))

  out <- leaflet::leaflet()
  out <- leaflet::addTiles(out)
  out <- leaflet::addCircleMarkers(out, lng = x[[lon_col]], lat = x[[lat_col]], popup = popup,
                                   radius = 8, stroke = FALSE, fillColor = "red", fillOpacity = 1,
                                   clusterOptions = leaflet::markerClusterOptions())
  return(out)
}





#' Barplots of the number of records for the most represented taxa
#'
#' Generate a multipanel plot where, for each taxonomic level,
#' a barplot represent the number of records available
#' in the reference database for the most represented taxa.
#'
#' @param x a reference database.
#' @param show_n an integer value indicating the number of taxa to show in each panel.
#'
#' @return A ggplot object.
#' @export
#'
refdb_plot_tax_barplot <- function(x, show_n = 10) {

  check_fields(x, c("taxonomy"))
  tax_cols <- attributes(x)$refdb$taxonomy

  dat <- tidyr::pivot_longer(x[tax_cols], tidyr::everything())
  dat <- dplyr::group_by(dat, .data$name, .data$value)
  dat <- dplyr::count(dat)
  dat <- dplyr::group_by(dat, .data$name)

  dat_na <- dplyr::filter(dat, is.na(.data$value))

  dat_top <- dplyr::filter(dat, !is.na(.data$value))
  dat_top <- dplyr::slice_max(dat_top, order_by = .data$n, n = show_n)
  dat_top <- dplyr::mutate(dat_top, val_ord = factor(.data$value,
                                                     levels = .data$value[order(.data$n, decreasing = TRUE)]))

  p <-
    ggplot2::ggplot(dat_top) +
    ggplot2::geom_col(ggplot2::aes(.data$val_ord, .data$n)) +
    ggplot2::facet_wrap(ggplot2::vars(.data$name), scales = "free") +
    ggplot2::ylab("Number of records") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1),
                   axis.title.x = ggplot2::element_blank())

  if(nrow(dat_na) > 0) {
    p <- p +
    ggplot2::geom_text(ggplot2::aes(Inf, Inf, label = paste("NA =", .data$n), hjust = 1.1, vjust = 1.5), size = 3, data = dat_na)
  }
  p
}
