
#' Compute the explained variance
#'
#' This function uses the standard deviation and the total variance
#' computed by Seurat when creating a PCA reduction dimension,
#' and computes the explained variance for each of the PCs.
#' It then adds the variance explained to the Seurat object,
#' under Seurat@@reductions$pca@@misc$variance_explained
#'
#' @param seur The Seurat object
#' @param reduction_name The name of the reduction to use. Default to "pca"
#'
#' @return The Seurat object with variance_explained in misc of dimreduc used
#' @export
varianceExplained <- function(seur, reduction_name = "pca") {
  pca <- seur@reductions[[reduction_name]]
  total_var <- pca@misc$total.variance
  eigValues <- (pca@stdev)^2
  varExplained <- eigValues/total_var
  seur@reductions[[reduction_name]]@misc$variance_explained <- varExplained
  return(seur)
}
