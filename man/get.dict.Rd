% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get.dict.R
\name{get.dict}
\alias{get.dict}
\title{Get the Ensembl dictionary}
\usage{

  get.dict(url = "https://www.ensembl.org/info/genome/stable_ids/prefixes.html",
  withdatasets = TRUE, host = "www.ensembl.org")
}
\arguments{
\item{url}{A character containing the URL of the official Ensembl site containing the prefix table. Defaults to the European mirror (www.ensembl.org/info/genome/stable_ids/prefixes.html)}

\item{withdatasets}{A logical indicating whether it should query Ensembl for dataset names. This is necessary for \code{get.orthoIDs} and \code{get.seq} functions. Defaults to TRUE.}

\item{host}{Ensembl host. Defaults to the European mirror.}
}
\value{
A data frame containing the Ensembl dictionary.
}
\description{
This function retrieves our creates an Ensembl dictionary, which includes species names in different formats, as well as IDs prefixes and dataset names for all available species names.
}
\examples{
get.dict()

}
