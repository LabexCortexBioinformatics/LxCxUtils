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

#' Count the number of cells
#'
#' Count the number of cells respecting all the conditions, and print a recap with the number of cells respecting each condition.
#'
#' @param seur Seurat object to use.
#' @param features Features to use.
#' @param assay Assay to use if features are genes. Default to "RNA".
#' @param slot Slot to use if features are genes. Default to "data".
#' @param value_to_compare Which value (or values) to compare the features to.
#' Must be either a single element or a list of the same length as features.
#' Default to 0.
#' @param comparison_type one of "inferior" or 1, "inferior or equal" or 2,
#' "equal" or 3, "superior or equal" or 4, "superior" or 5.
#' Can be a list of them, the same length as features. Default to "superior".
#' @param verbose Whether or not to print the summary dataframe. Default to TRUE.
#'
#' @return The number of cells respecting all the conditions.
#' @export
WhichCells2 <- function(seur,
                         features,
                         assay = "RNA",
                         slot = "data",
                         value_to_compare = 0,
                         comparison_type = "superior",
                         verbose = TRUE){

  comparison = unlist(lapply(comparison_type, function(x)
                             switch(x,
                                    "inferior" = `<`,
                                    "inferior or equal" = `<=`,
                                    "equal" = `==`,
                                    "superior or equal" = `>=`,
                                    "superior" = `>`,
                                    "not equal" = `!=`)))

  if (length(comparison_type) == 1){
    comparison = rep(comparison, times = length(features))
  }
  names(comparison) = features

  if (length(value_to_compare) == 1){
    value_to_compare <- rep(value_to_compare, times = length(features))
  }
  names(value_to_compare) <- features

  def <- DefaultAssay(seur)
  DefaultAssay(seur) <- assay
  data <- FetchData(seur, c(features, "ident"), slot = slot)
  DefaultAssay(seur) <- def
  df <- data.frame(matrix(ncol = 0, nrow = 1))
  row.names(df) <- "number of cells"

  cells <- lapply(features, function(f){
    name <- grep(f, colnames(data), value = TRUE)[1]
    comp <- comparison[[f]]
    in_ <- comp(data[[name]], value_to_compare[f])
    row.names(data[in_,])
  })
  if (verbose){
    for (i in 1:length(cells)){
      df[[colnames(data)[i]]] <- length(unlist(cells[i]))
    }
    print(df)
  }
  if (length(features) > 1){
    Reduce(intersect, cells)
  } else {
    unlist(cells)
  }
}


#' Count the number of cells respecting the conditions
#'
#' Apply WhichCells2 and return the length.
#'
#' @param seur Seurat object to use.
#' @param features Features to use.
#' @param assay Assay to use if features are genes. Default to "RNA".
#' @param slot Slot to use if features are genes. Default to "data".
#' @param value_to_compare Which value (or values) to compare the features to.
#' Must be either a single element or a list of the same length as features.
#' Default to 0.
#' @param comparison_type one of "inferior" or 1, "inferior or equal" or 2,
#' "equal" or 3, "superior or equal" or 4, "superior" or 5.
#' Can be a list of them, the same length as features. Default to "superior".
#' @param verbose Whether or not to print the summary dataframe. Default to TRUE.
#'
#' @return The number of cells respecting all the conditions.
#' @export
howManyCells <- function(seur,
                         features,
                         assay = "RNA",
                         slot = "data",
                         value_to_compare = 0,
                         comparison_type = "superior",
                         verbose = TRUE) {

  length(WhichCells2(seur,
                     features,
                     assay,
                     slot,
                     value_to_compare,
                     comparison_type,
                     verbose)
         )
}

#' Set colors
#'
#' @param cols NULL or a vector of colors (named or not)
#' @param group.by column with duplicated names of colors
#'
#' @return colors named
setCols <- function(cols, group.by){
  if (is.factor(group.by)){
    unique_val <- levels(group.by)
  } else {
    unique_val <- unique(group.by)
  }
  if (is.null(cols)) {
    cols <- scales::hue_pal()(length(unique_val))
  }
  if (is.null(names(cols))) {
    names(cols) <- unique_val
  }
  return(cols)
}



