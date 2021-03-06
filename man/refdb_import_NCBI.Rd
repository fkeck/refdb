% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import_NCBI.R
\name{refdb_import_NCBI}
\alias{refdb_import_NCBI}
\title{Download and import NCBI Nucleotide records}
\usage{
refdb_import_NCBI(query, full = FALSE)
}
\arguments{
\item{query}{a character string with the query}

\item{full}{a logical. If FALSE (the default), only a subset of
the most important fields is included in the result.}
}
\value{
A tibble.
}
\description{
This function allows to search and download data
from the the NCBI Nucleotide database.
Additionally it uses the NCBI Taxonomy database
to get the sequence taxonomic classification.
}
\details{
This function uses several functions of the \pkg{rentrez}
package to interface with the NCBI's EUtils API.
}