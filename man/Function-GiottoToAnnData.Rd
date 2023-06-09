% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GiottoToAnndata.R
\name{Function GiottoToAnnData}
\alias{Function GiottoToAnnData}
\alias{GiottoToAnnData}
\title{GiottoToAnnData}
\usage{
GiottoToAnnData(
  object = NULL,
  outpath = NULL,
  assays = NULL,
  reductions = NULL,
  images = NULL,
  markersDF = NULL
)
}
\arguments{
\item{object}{giotto object used to make anndata}

\item{outpath}{anndata output path}

\item{assays}{name of assays to export}

\item{reductions}{a vector of reduction names used to export}

\item{images}{name of image used to export}

\item{markersDF}{a named list of marker df vars}
}
\description{
Coveret Giotto To Anndata
}
