% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seurat_plot.R
\name{seuratBoxPlot}
\alias{seuratBoxPlot}
\title{\subsection{Boxplot Function}{
}}
\usage{
seuratBoxPlot(
  seur,
  features,
  group.by = "ident",
  assay = "RNA",
  slot = "data",
  ncol = 2,
  cols = NULL,
  pt.size = 1,
  pt.alpha = 1,
  remove_legend = FALSE,
  remove_x_text = FALSE
)
}
\arguments{
\item{seur}{The Seurat object within which the data will be fetched.}

\item{features}{Vector of features to use.
Can be genes, or any other continuous column accessible by Seurat's FetchData.}

\item{group.by}{Metadata column used to make the groups in the boxplot.
Default to "ident"}

\item{assay}{Assay used to Fetch the data from. Will replace the Default Assay.
Default to "RNA"}

\item{slot}{Slot used within the assay the fetch data from. Default to "data"}

\item{ncol}{Number of column in the gridExtra object. Default to 2}

\item{pt.size}{Size of the points to plot.
Either set this to 0 or pt.alpha to 0 to remove points. Default to 1}

\item{pt.alpha}{Transparency of the points.
Either set this to 0 or pt.size to 0 to remove points. Default to 1}

\item{remove_legend}{Whether to remove the color legend. Default to FALSE}

\item{remove_x_text}{Whether to remove x-axis legend. Default to FALSE}
}
\value{
A gridExtra object containing as much plots as there are elements in the \code{features} parameter.
If there is only one feature, it will return the ggplot object.
}
\description{
This function aims at replacing the VlnPlot function from Seurat with box plots.
}
