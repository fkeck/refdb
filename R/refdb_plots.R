

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
  res_tidy <- dplyr::group_by(res_tidy, taxonomy_rank)
  res_tidy <- dplyr::mutate(res_tidy, freq_by_rank = n/sum(n))
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
#' @param x
#' @param cols a vector of column names to include in the treemap.
#' If not provided, the function tries to find a relevant subset of columns.
#' @param freq_labels Use these parameters to adjust the number of labels printed
#'
#' @details
#' TODO
#'
#' Only the columns provided in the cols argument are represented
#' in the treemap.
#' Large labels are printed for the highest rank, while light
#' text labels are printed for the lowest rank
#' Intermediate ranks are drawn but their names are not shown.
#'
#' Use internal function to make a graph and design your plot
#'
#' @return
#' A \pkg{ggplot2} (\pkg{ggraph}) object. This means the plot can be further customized
#' using ggplot2 compatible functions.
#' @export
#'
#' @examples
refdb_plot_treemap <- function(x, cols = NULL, freq_labels = c(0.01, 0.003)) {

  check_fields(x, what = "taxonomy")
  col_tax <- attributes(x)$refdb_fields$taxonomy


  if(is.null(cols)) {
    cols <- apply(x[, col_tax], 2, function(x) length(unique(x)))
    cols <- names(cols[cols > 50 & cols < 5000])
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

  ggraph::ggraph(g, 'treemap', weight =  n) +
    ggraph::geom_node_tile(ggplot2::aes(fill = taxo_nth_level(name, 1),
                                        alpha = n),
                           size = 0, show.legend = FALSE) +
    ggraph::geom_node_tile(ggplot2::aes(size = taxonomy_depth),
                           colour = 'white') +
    ggraph::geom_node_text(ggplot2::aes(label = ifelse(taxonomy_depth == max(taxonomy_depth) & freq_by_rank > freq_labels[2], terminal, NA_character_)),
                           size = 3, color = "white", alpha = 0.8) +
    ggraph::geom_node_label(ggplot2::aes(label = ifelse(taxonomy_depth == 1 & freq_by_rank > freq_labels[1], terminal, NA_character_))) +
    ggplot2::scale_alpha(range = c(1, 0.5), guide = 'none') +
    ggplot2::scale_size(range = c(4, 0.2), guide = 'none')

}


#
#
#
# ggraph::ggraph(g, layout = "igraph", algorithm = "fr") +
#   ggraph::geom_edge_link() +
#   ggraph::geom_node_point(ggplot2::aes(size = n, color = as.factor(taxonomy_depth)))
#
# ggraph::ggraph(g, layout = "circlepack", circular = TRUE) +
#   ggraph::geom_edge_link() +
#   ggraph::geom_node_point(ggplot2::aes(size = n, color = as.factor(taxonomy_depth))) +
#   ggplot2::coord_fixed()
#
#
# ggraph::ggraph(g, layout = "unrooted", daylight = T) +
#   ggraph::geom_edge_link() +
#   ggraph::geom_node_point(ggplot2::aes(size = refdb_n, color = as.factor(refdb_depth))) +
#   ggplot2::coord_fixed()
#
#
# ggraph::ggraph(g, layout = "auto") +
#   ggraph::geom_edge_link() +
#   ggraph::geom_node_point()
#

#
#
#
# ggraph::ggraph(g, 'dendrogram', circular = TRUE)) +
#   ggraph::geom_edge_elbow() +
#   ggraph::geom_node_point(ggplot2::aes(size = refdb_n, color = as.factor(refdb_depth))) +
#   ggplot2::coord_fixed()
#
#
# ggraph::ggraph(g, 'dendrogram') +
#   ggraph::geom_edge_elbow()
#
# ggraph::ggraph(g, 'partition') +
#   ggraph::geom_node_tile(size = 0.25)
#
#
