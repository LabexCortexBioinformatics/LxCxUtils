% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seurat_plot.R
\name{QC_barplot}
\alias{QC_barplot}
\title{Plot features as barplot}
\usage{
QC_barplot(
  seur,
  features,
  group.by = "orig.ident",
  cols = NULL,
  wrap = FALSE,
  assay = "RNA",
  slot = "data",
  bins = 100,
  ncol = 2
)
}
\arguments{
\item{seur}{Seurat object to fetch data from.}

\item{features}{Features to use. Use nFeature_RNA, nCount_RNA, percent.mt or percent.ribo for QCs.}

\item{group.by}{Which variable to group by. Default to "orig.ident"}

\item{cols}{Which colors to use. NULL will use hue_pal. Default to NULL}

\item{wrap}{Whether or not to use the facet_wrap from ggplot on the group.by groups. Default to FALSE}

\item{assay}{Which assay to get data from, if features are genes. Default to "RNA"}

\item{slot}{Which slot in the assay to get data from, if features are genes. Default to "data"}

\item{bins}{Number of bins to use. A high number will have better precision. Default to 100}

\item{ncol}{Number of column if several features. Default to 2}
}
\value{
A gridExtra object if several features, or a ggplot object if only one feature.
}
\description{
Does the same as Violin or Box plots, but as an histogram.
}
