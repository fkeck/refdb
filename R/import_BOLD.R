

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

  if(query_stats$order$count > 1 & query_stats$total_records > 5000) {

    # split by order strategy

  } else {
    cat("Downloading", query_stats$total_records, "sequences from BOLD...\n")

    recs <- bold::bold_seqspec(taxon = taxon,
                               ids = ids,
                               bin = bin,
                               container = container,
                               institutions = institutions,
                               researchers = researchers,
                               geo = geo)

    out <- tibble::as_tibble(recs)
  }

  out <- tibble::tibble(source = "BOLD", out)

  if (full == FALSE) {
    out <- dplyr::select(out, source, sequenceID, markercode,
                         phylum_name, class_name, order_name,
                         family_name, subfamily_name, genus_name,
                         species_name, subspecies_name, # Correspond to NCBI organism
                         nucleotides,
                         country, province_state, lat, lon)
  }

  # Empty strings as NA
  out[out == ""] <- NA

  #
  if (ncbi_taxo == TRUE) {

    bt_taxo_cols <- c("genus_name", "family_name", "order_name", "class_name", "phylum_name")
    ncbi_taxo_cols <- c("genus", "family", "order", "class", "phylum")

    bt <- out[, bt_taxo_cols]
    bt <- dplyr::distinct(bt)

    ncbi_ids <- vector("character", nrow(bt))
    for(i in seq_along(ncbi_ids)) {

      x <- unlist(bt[i, ])
      ncbi_ids[i] <- NA
      lvl = 1

      while (is.na(ncbi_ids[i])) {

        if(is.na(x[lvl])) {
          ncbi_ids[i] <- NA
        } else {
          cat("\rProcessing:", x[lvl], rep(" ", 40))
          suppressMessages(
            uid <- taxize::get_uid_(x[lvl], messages = FALSE)[[1]]$uid
          )

          # Deal with no match
          if(length(uid) == 0) {
            ncbi_ids[i] <- NA
          }

          # Deal with multiple matches
          if(length(uid > 1)) {
            tax <- get_ncbi_taxonomy(uid)
            score <- apply(tax[ncbi_taxo_cols], 1,
                           function(y) sum(y == x, na.rm = TRUE)
            )
            ncbi_ids[i] <- tax$id[which.max(score)]
          }
        }
        lvl = lvl + 1
      }
    }

    bt_taxo <- get_ncbi_taxonomy(ncbi_ids)
    bt_taxo <- dplyr::bind_cols(bt, bt_taxo)
    out <- dplyr::left_join(out, bt_taxo, by = bt_taxo_cols)

  }

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

  return(out)
}


