% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/refdb_ncbi_taxonomy.R
\name{refdb_set_ncbitax}
\alias{refdb_set_ncbitax}
\title{Replace the current taxonomy using the NCBI Taxonomy database}
\usage{
refdb_set_ncbitax(x)
}
\arguments{
\item{x}{a reference database (tibble) with one or several columns
giving the taxonomy of each record and explicitly
indicated in the field taxonomy. See \link{refdb_set_fields}.}
}
\value{
The reference database with the NCBI taxonomy
for the genus level and higher ranks.
(the original taxonomy above the genus level is removed).
}
\description{
Replace the current taxonomy using the NCBI Taxonomy database
}
