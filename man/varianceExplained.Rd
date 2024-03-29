% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{varianceExplained}
\alias{varianceExplained}
\title{Compute the explained variance}
\usage{
varianceExplained(seur, reduction_name = "pca")
}
\arguments{
\item{seur}{The Seurat object}

\item{reduction_name}{The name of the reduction to use. Default to "pca"}
}
\value{
The Seurat object with variance_explained in misc of dimreduc used
}
\description{
This function uses the standard deviation and the total variance
computed by Seurat when creating a PCA reduction dimension,
and computes the explained variance for each of the PCs.
It then adds the variance explained to the Seurat object,
under Seurat@reductions$pca@misc$variance_explained
}
