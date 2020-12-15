

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

  if (ncbi_taxo) {
    out <- refdb_ncbi_taxonomy(out)
  }

  return(out)
}


