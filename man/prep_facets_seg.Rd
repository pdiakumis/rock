% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cnv_facets.R
\name{prep_facets_seg}
\alias{prep_facets_seg}
\title{Read FACETS CNV Segments}
\usage{
prep_facets_seg(facets)
}
\arguments{
\item{facets}{Path to FACETS \code{emcncf} text file.}
}
\value{
A dataframe (\code{tibble}) with the following columns:
\itemize{
\item chrom: chromosome
\item start: start coordinate
\item end: end coordinate
\item tot_cn: total copy number estimate
}
}
\description{
Reads the output by the FACETS \code{emcncf} function and exports the \code{cncf}
CNV segment coordinates.
}
\examples{
cn <- system.file("extdata", "HCC2218_facets_cncf.tsv", package = "pebbles")
prep_facets_seg(cn)

}
