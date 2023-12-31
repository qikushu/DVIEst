% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/01_functions.R
\name{getEstParam}
\alias{getEstParam}
\title{Get estimated parameters of a DVI model}
\usage{
getEstParam(object)
}
\arguments{
\item{object}{a DviEstl object created by \code{\link[=estDVIparam]{estDVIparam()}}.}
}
\value{
a numeric vector of estimated alpha, beta, and g.
}
\description{
Extracted esetimted parameter values from the DviEst object.
}
