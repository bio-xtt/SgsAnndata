% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SeuratToAnndata2.R
\name{Function SeuratToAnndata}
\alias{Function SeuratToAnndata}
\alias{SeuratToAnndata}
\title{SeuratToAnndata}
\usage{
SeuratToAnndata(
  object,
  outpath,
  main.assays,
  assays = NULL,
  groups = NULL,
  reductions = NULL,
  images = NULL,
  markersDF = NULL
)
}
\arguments{
\item{object}{seurat object}

\item{outpath}{outputpath of the anndata object}

\item{main.assays}{name of assay used to build adata.X}

\item{assays}{a vector of assay names to export in adata}

\item{groups}{vector of groups to export,if null export all}

\item{reductions}{vector of reduction names to export,if null export all}

\item{images}{vector of images names to export,if null export all}

\item{markersDF}{a dataframe of the marker files,need the colnames in c(gene, cluster, ...) order}
}
\description{
SeuratToAnndata
}
