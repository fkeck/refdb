

#' Download and import NCBI Nucleotide records
#'
#' This function allows to search and download data
#' from the the NCBI Nucleotide database.
#' Additionally it uses the NCBI Taxonomy database
#' to get the sequence taxonomic classification.
#'
#' @param query a character string with the query
#' @param full a logical. If FALSE (the default), only a subset of
#' the most important fields is included in the result.
#'
#' @details
#' This function uses several functions of the \pkg{rentrez}
#' package to interface with the NCBI's EUtils API.
#'
#' @return A tibble.
#' @export
#'
refdb_import_NCBI <- function(query, full = FALSE) {

  ff <- tempfile("refdb_NCBI_", fileext = ".csv")

  req <- rentrez::entrez_search(db = "nuccore",
                                term = query,
                                use_history = TRUE)

  cat("Downloading", req$count, "sequences from NCBI...\n")

  pb <- utils::txtProgressBar(0, req$count, 0, style = 3)

  # Main loop to download and write the data
  for(seq_start in seq(0, req$count, 200)){

    recs <- rentrez::entrez_fetch(db = "nuccore",
                                  web_history = req$web_history,
                                  rettype = "gb", retmode = "xml",
                                  retmax = 200, retstart = seq_start)

    NCBI_xml <- xml2::read_xml(recs)
    NCBI_xml <- xml2::xml_children(NCBI_xml)

    NCBI_table <- make_ncbi_table(NCBI_xml)

    taxo_id <- xml2::xml_text(xml2::xml_find_all(NCBI_xml, './/GBQualifier_name[text()="db_xref"]/following-sibling::GBQualifier_value'))
    taxo_id <- taxo_id[stringr::str_detect(taxo_id, "taxon:[0-9]+")]
    taxo_id <- stringr::str_extract(taxo_id, "(?<=taxon:)[0-9]+")
    taxo_id <- tibble::tibble(taxonomy = NCBI_table$taxonomy, id = taxo_id)
    taxo_id <- taxo_id[!duplicated(taxo_id$taxonomy), ]
    taxo_id <- dplyr::left_join(taxo_id,
                                get_ncbi_taxonomy(taxo_id$id),
                                by = "id")

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

    utils::setTxtProgressBar(pb, seq_start)

  }

  utils::setTxtProgressBar(pb, req$count)

  out <- readr::read_csv(ff, col_types = readr::cols())

  out <- refdb_set_fields(out,
                          source = "source",
                          id = "id",
                          taxonomy = c(superkingdom = "superkingdom",
                                       kingdom = "kingdom",
                                       phylum = "phylum",
                                       subphylum = "subphylum",
                                       class = "class",
                                       subclass = "subclass",
                                       infraclass = "infraclass",
                                       order = "order",
                                       suborder = "suborder",
                                       infraorder = "infraorder",
                                       superfamily = "superfamily",
                                       family = "family",
                                       genus = "genus",
                                       species = "species"),
                          sequence = "sequence",
                          marker = "gene")

  file.remove(ff)
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

  taxo <- rentrez::entrez_fetch("taxonomy", id = id, rettype = "xml")
  taxo_xml <- xml2::read_xml(taxo)
  taxo_xml <- xml2::xml_children(taxo_xml)

  taxo_table <- tibble::tibble(
    id = id,
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
    genus = xml_extract(taxo_xml, './/Rank[text()="genus"]/preceding-sibling::ScientificName')
  )

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
