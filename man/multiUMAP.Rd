% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seurat_plot.R
\name{multiUMAP}
\alias{multiUMAP}
\title{Compute multiple UMAPs}
\usage{
multiUMAP(
  seur,
  dims,
  out = TRUE,
  save = FALSE,
  add_to_object = FALSE,
  group.by = "ident",
  ncol = 3,
  pt.size = 0.5,
  file = paste0(getwd(), "/multiUMAP"),
  fig_size = "auto",
  dpi = 200,
  save_as = "png",
  verbose = TRUE
)
}
\arguments{
\item{seur}{Seurat object to use.}

\item{dims}{Dimensions to use. Must be an integer vector or list.
Can be cast with a range, for example 5:10 will compute UMAPs 5, 6, 7, 8, 9, 10.}

\item{out}{Whether to return the aggregated plot. Default to TRUE}

\item{save}{Whether to save the aggregated plot. Default to FALSE}

\item{add_to_object}{Whether to add every UMAP computed in the object.
If TRUE, the UMAPs will be saved under "UMAP_dim" followed by the dimnesion
(for exemple "UMAP_dim5" for dimension 5). If FALSE, the UMAPs will be
computed under "temp" and then suppressed. Default to FALSE}

\item{group.by}{Which column of the metadata to group cells by. Default to ident, which will take the active ident.}

\item{ncol}{Number of columns in the aggregated plot. Default to 3}

\item{file}{The file name and path if image is saved. Default to working_directory/multiUMAP.png}

\item{fig_size}{The size of the saved image. Either "auto", which will try to use the best image size automatically,
or a vector of int, the first one being the width in pixel and second one the height in pixel.
Default is "auto"}

\item{dpi}{The dpi of the saved image. Must be an int. Default to 200}

\item{save_as}{The type of file to save. Default to "png"}

\item{verbose}{Default to TRUE}

\item{pt_size}{Size of points to use. Default to 0.5}
}
\value{
a gridExtra object if TRUE, else nothing
}
\description{
Will compute multiple umaps, and aggregate their dimplots in one plot.
This plot can then be returned or saved as image (or both).
}
