% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/refdb_clean.R
\name{refdb_clean_tax_NA}
\alias{refdb_clean_tax_NA}
\title{Convert missing taxonomic names to NA}
\usage{
refdb_clean_tax_NA(x, cols = NULL)
}
\arguments{
\item{x}{a reference database.}

\item{cols}{an optional vector of column names.
If \code{NULL} (default), the function is applied to the columns
associated with the \code{taxonomy} and \code{organism} fields.}
}
\value{
A reference database.
}
\description{
Convert missing taxonomic names to NA
}
