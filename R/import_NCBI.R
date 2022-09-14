

#' Download and import NCBI Nucleotide records
#'
#' This function allows to search and download data
#' from the the NCBI Nucleotide database.
#' Additionally it uses the NCBI Taxonomy database
#' to get the sequence taxonomic classification.
#'
#' @param query a character string with the query.
#' @param full a logical. If FALSE (the default), only a subset of
#' the most important fields is included in the result.
#' @param max_seq_length a numeric giving the maximum length of sequences
#' to retrieve. Useful to exclude complete genomes.
#' @param seq_bin number of sequences to download at once.
#'
#' @details
#' This function uses several functions of the \pkg{rentrez}
#' package to interface with the NCBI's EUtils API.*
#'
#' @section Errors:
#' \code{Error in curl::curl_fetch_memory(url, handle = handle) :
#' transfer closed with outstanding read data remaining} \cr
#' This error seems to appear with long sequences.
#' You can try to decrease \code{max_seq_length} to exclude them.
#'
#' @return A tibble.
#'
#' @examples
#' \dontrun{
#' silo_ncbi <- refdb_import_NCBI("Silo COI")
#' }
#'
#' @export
#'
refdb_import_NCBI <- function(query,
                              full = FALSE,
                              max_seq_length = 10000,
                              seq_bin = 200) {

  ff <- tempfile("refdb_NCBI_", fileext = ".csv")
  fx <- tempfile("refdb_NCBI_", fileext = ".xml")

  query <- paste0(query, ' AND ( "0"[SLEN] : "', max_seq_length, '"[SLEN] )')

  req <- rentrez::entrez_search(db = "nuccore",
                                term = query,
                                use_history = TRUE)

  if(req$count == 0) {
    cat("No sequence found\n")
    return(NULL)
  }

  cat("Downloading", req$count, "sequences from NCBI...\n")

  # Main loop to download and write the data
  for(seq_start in seq(0, req$count, seq_bin)){

    recs <- entrez_fetch_retry(db = "nuccore",
                               web_history = req$web_history,
                               rettype = "gb", retmode = "xml",
                               retmax = seq_bin, retstart = seq_start,
                               delay_retry = 60, n_retry = 50)


    if(is.na(recs)) {
      next
      # There is a risk of loss of data here (we drop seq_bin records)
      # Need to track these cases which seems to be linked
      # to NCBI empty records included in the search results
    }

    readr::write_lines(recs, file = fx, append = FALSE)
    NCBI_xml <- xml2::read_xml(fx)
    NCBI_xml <- xml2::xml_children(NCBI_xml)

    NCBI_table <- make_ncbi_table(NCBI_xml)

    taxo_id <- xml2::xml_text(xml2::xml_find_all(NCBI_xml, './/GBQualifier_name[text()="db_xref"]/following-sibling::GBQualifier_value'))
    taxo_id <- taxo_id[stringr::str_detect(taxo_id, "taxon:[0-9]+")]
    taxo_id <- stringr::str_extract(taxo_id, "(?<=taxon:)[0-9]+")
    taxo_id <- tibble::tibble(taxonomy = NCBI_table$taxonomy, id = taxo_id)
    taxo_id <- taxo_id[!duplicated(taxo_id$taxonomy), ]
    gtax <- get_ncbi_taxonomy_retry(taxo_id$id, delay_retry = 60, n_retry = 50)
    taxo_id <- dplyr::left_join(taxo_id, gtax[, -ncol(gtax)], by = "id")

    NCBI_table <- dplyr::left_join(NCBI_table, taxo_id,
                                   by = "taxonomy",
                                   suffix = c("", "_taxonomy"))

    NCBI_table <- tibble::tibble(source = "NCBI", NCBI_table)
    NCBI_table <- dplyr::mutate(NCBI_table, species = .data$organism)
    # Maybe clean species column here

    if(full == FALSE) {
      NCBI_table <- dplyr::select(NCBI_table,
                                  .data$source,
                                  .data$id,
                                  .data$gene,
                                  .data$sequence,
                                  .data$superkingdom,
                                  .data$kingdom,
                                  .data$phylum,
                                  .data$subphylum,
                                  .data$class,
                                  .data$subclass,
                                  .data$infraclass,
                                  .data$order,
                                  .data$suborder,
                                  .data$infraorder,
                                  .data$superfamily,
                                  .data$family,
                                  .data$genus,
                                  .data$species,
                                  .data$country_location,
                                  .data$lat_lon)
    }


    # Write header
    if(seq_start == 0) {
      readr::write_csv(NCBI_table[0, ], file = ff)
    }

    readr::write_csv(NCBI_table,
                     file = ff,
                     append = TRUE,
                     col_names = FALSE)

    cat("\r > ", seq_start + nrow(NCBI_table),
        " (",
        round((seq_start + nrow(NCBI_table))/req$count * 100, digits = 1),
        "%) ",
        "sequences downloaded.", sep = "")

  }


  out <- readr::read_csv(ff, col_types = readr::cols())
  # Process geographic coordinates
  out <- process_geo_ncbi(out)

  # Set fields
  out <- refdb_set_fields_NCBI(out)
  out <- refdb_set_fields(out,
                          latitude = "latitude",
                          longitude = "longitude")

  file.remove(ff, fx)
  return(out)
}




#' Get NCBI taxonomy
#'
#' Download and parse NCBI taxonomy records
#'
#' @param id A vector of id for records in the NCBI Taxonomy database.
#'
#' @return A tibble with each row corresponding to an id and each column
#' to a taxonomic level.
get_ncbi_taxonomy <- function(id) {

  ids <- split(id, ceiling(seq_along(id)/100))

  taxo_table <- lapply(ids, function(x) {

    taxo <- entrez_fetch_retry("taxonomy", id = x, rettype = "xml")

    taxo_xml <- xml2::read_xml(taxo)
    taxo_xml <- xml2::xml_children(taxo_xml)

    taxo_table <- tibble::tibble(
      id = x,
      superkingdom = xml_extract(taxo_xml, './/Rank[text()="superkingdom"]/preceding-sibling::ScientificName'),
      kingdom = xml_extract(taxo_xml, './/Rank[text()="kingdom"]/preceding-sibling::ScientificName'),
      phylum = xml_extract(taxo_xml, './/Rank[text()="phylum"]/preceding-sibling::ScientificName'),
      subphylum = xml_extract(taxo_xml, './/Rank[text()="subphylum"]/preceding-sibling::ScientificName'),
      class = xml_extract(taxo_xml, './/Rank[text()="class"]/preceding-sibling::ScientificName'),
      subclass = xml_extract(taxo_xml, './/Rank[text()="subclass"]/preceding-sibling::ScientificName'),
      infraclass = xml_extract(taxo_xml, './/Rank[text()="infraclass"]/preceding-sibling::ScientificName'),
      order = xml_extract(taxo_xml, './/Rank[text()="order"]/preceding-sibling::ScientificName'),
      suborder = xml_extract(taxo_xml, './/Rank[text()="suborder"]/preceding-sibling::ScientificName'),
      infraorder = xml_extract(taxo_xml, './/Rank[text()="infraorder"]/preceding-sibling::ScientificName'),
      superfamily = xml_extract(taxo_xml, './/Rank[text()="superfamily"]/preceding-sibling::ScientificName'),
      family = xml_extract(taxo_xml, './/Rank[text()="family"]/preceding-sibling::ScientificName'),
      genus = xml_extract(taxo_xml, './/Rank[text()="genus"]/preceding-sibling::ScientificName'),
      species = xml_extract(taxo_xml, './/Rank[text()="species"]/preceding-sibling::ScientificName')
    )

    return(taxo_table)
  })

  taxo_table <- dplyr::bind_rows(taxo_table)

  return(taxo_table)
}




#' Parse NCBI XML and make a table
#'
#' @param x A XML nodeset.
#'
#' @return A tibble.
#'
make_ncbi_table <- function(x) { # Add all possible fields Gene is missing!
  tibble::tibble(
    id = xml_extract(x, "./GBSeq_primary-accession"),
    organism = xml_extract(x, "./GBSeq_organism"),
    taxonomy = xml_extract(x, "./GBSeq_taxonomy"),
    sequence = xml_extract(x, "./GBSeq_sequence"),
    create_date = xml_extract(x, "./GBSeq_create-date"),
    update_date = xml_extract(x, "./GBSeq_update-date"),
    definition = xml_extract(x, "./GBSeq_definition"),
    other_seqids = xml_extract(x, "./GBSeq_other-seqids"),
    moltype = xml_extract(x, "./GBSeq_moltype"),
    length = xml_extract(x, "./GBSeq_length"),
    topology = xml_extract(x, "./GBSeq_topology"),
    gene = xml_extract(x, './/GBQualifier_name[text()="gene"]/following-sibling::GBQualifier_value'),
    product = xml_extract(x, './/GBQualifier_name[text()="product"]/following-sibling::GBQualifier_value'),
    country_location = xml_extract(x, './/GBQualifier_name[text()="country"]/following-sibling::GBQualifier_value'),
    lat_lon = xml_extract(x, './/GBQualifier_name[text()="lat_lon"]/following-sibling::GBQualifier_value'),
    isolation_source = xml_extract(x, './/GBQualifier_name[text()="isolation_source"]/following-sibling::GBQualifier_value'),
    specimen_voucher = xml_extract(x, './/GBQualifier_name[text()="specimen_voucher"]/following-sibling::GBQualifier_value'),
    reference_authors = sapply(xml2::xml_find_first(x, "./GBSeq_references/GBReference/GBReference_authors"),
                                    function(x) paste(xml2::xml_text(xml2::xml_find_all(x, "./GBAuthor")), collapse = ", ")),
    reference_title = xml_extract(x, "./GBSeq_references/GBReference/GBReference_title"),
    reference_journal = xml_extract(x, "./GBSeq_references/GBReference/GBReference_journal")
  )
}


#' Extract XML elements
#'
#' Combine xml_find_first and xml_text to extract elements.
#'
#' @param x A document, node, or node set.
#' @param xpath A string containing a xpath expression.
#'
#' @return A character vector, the same length as x.
#'
xml_extract <- function(x, xpath) {
  xml2::xml_text(xml2::xml_find_first(x, xpath))
}



#' Process coordinate column returned by NCBI
#'
#' @param x NCBI dataframe.
#' @param col column name containing geographical coordinates.
#'
#' @return NCBI dataframe.
#'
process_geo_ncbi <- function(x, col = "lat_lon") {
  x <- tidyr::separate(x, col = col, into = c("latitude", "longitude"), sep = "(?<=[A-Z]) ")
  lat_sign <- ifelse(stringr::str_detect(x$latitude, "N"), 1, -1)
  x$latitude <- as.numeric(stringr::str_remove(x$latitude, " [A-Z]"))
  x$latitude <- x$latitude * lat_sign

  lon_sign <- ifelse(stringr::str_detect(x$longitude, "E"), 1, -1)
  x$longitude <- as.numeric(stringr::str_remove(x$longitude, " [A-Z]"))
  x$longitude <- x$longitude * lon_sign

  return(x)
}




# Retry entrez_fetch
# Will retry every delay_retry seconds for n_retry times
entrez_fetch_retry <- function(..., delay_retry = 60, n_retry = 20) {

  res <- "error"

  while (res == "error" & n_retry > 0) {

    res <- tryCatch({
      Sys.sleep(0.1)
      # Abort if slower than 30 bytes/sec during 60 seconds
      # To avoid httr dependency we construct the config manually (not very safe)
      # Replace httr::config(low_speed_limit = 30L, low_speed_time = 60L)
      httr_conf <- structure(list(method = NULL, url = NULL, headers = NULL, fields = NULL,
                                  options = list(low_speed_limit = 30L, low_speed_time = 60L),
                                  auth_token = NULL, output = NULL), class = "request")


      rentrez::entrez_fetch(..., config = httr_conf)
    },
    error = function(cond) {
      message("\nSomething went wrong:")
      message(cond)
      message("\n")
      for (i in delay_retry:0) {cat("\rRetrying in", i, "s.  "); Sys.sleep(1)}
      cat("\n")
      return("error")
    }
    )
    n_retry <- n_retry - 1
  }

  if(res == "error") {
    stop("All attempts failed.")
  } else {
    return(res)
  }
}


# Another layer of security for large request
# and which proved to be useful
get_ncbi_taxonomy_retry <- function(id, delay_retry = 60, n_retry = 20) {

  res <- "error"

  while (identical(res, "error") & n_retry > 0) {

    res <- tryCatch({
      Sys.sleep(0.1)
      get_ncbi_taxonomy(id)
    },
    error = function(cond) {
      message("\nSomething went wrong:")
      message(cond)
      message("\n")
      for (i in delay_retry:0) {cat("\rRetrying in", i, "s.  "); Sys.sleep(1)}
      cat("\n")
      return("error")
    })
    n_retry <- n_retry - 1
  }

  if(identical(res, "error")) {
    stop("All attempts failed.")
  } else {
    return(res)
  }
}
