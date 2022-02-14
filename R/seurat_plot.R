library(Seurat)
library(ggplot2)


#' ggplot2 theme with white background and axis lines
#'
#' @return ggplot2 theme
theme_white_bg <- function(){
  theme(panel.background = element_blank(),
                 axis.line = element_line())
}

#' ## Boxplot Function
#' This function aims at replacing the VlnPlot function from Seurat with box plots.
#'
#' @param seur The Seurat object within which the data will be fetched.
#' @param features Vector of features to use.
#'     Can be genes, or any other continuous column accessible by Seurat's FetchData.
#' @param group.by Metadata column used to make the groups in the boxplot.
#'     Default to "ident"
#' @param assay Assay used to Fetch the data from. Will replace the Default Assay.
#'     Default to "RNA"
#' @param slot Slot used within the assay the fetch data from. Default to "data"
#' @param ncol Number of column in the gridExtra object. Default to 2
#' @param pt.size Size of the points to plot.
#'     Either set this to 0 or pt.alpha to 0 to remove points. Default to 1
#' @param pt.alpha Transparency of the points.
#'     Either set this to 0 or pt.size to 0 to remove points. Default to 1
#' @param remove_legend Whether to remove the color legend. Default to FALSE
#'
#' @return A gridExtra object containing as much plots as there are elements in the `features` parameter.
#'     If there is only one feature, it will return the ggplot object.
#' @export
seuratBoxPlot <- function(seur,
                    features,
                    group.by = "ident",
                    assay = "RNA",
                    slot = "data",
                    ncol = 2,
                    pt.size = 1,
                    pt.alpha = 1,
                    remove_legend = TRUE) {

  DefaultAssay(seur) <- assay
  plot_data <- FetchData(seur, vars = c(features, group.by), slot = slot)
  grobs <- list()

  for (f in features) {
    p <- ggplot(plot_data, aes(x = .data[[group.by]], y = .data[[f]], fill = .data[[group.by]])) +
      geom_boxplot(outlier.alpha = 0) +
      geom_jitter(size = ifelse(pt.size<=0, -1, pt.size), alpha = pt.alpha) +
      theme_white_bg() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, colour = "black"),
            plot.title = element_text(face="bold")) +
      labs(y = f, title = f)

    if (remove_legend) {
      p <- p + theme(legend.position = "none")
    }
    grobs[[f]] <- p
  }

  if (length(features > 1)) {
    gridExtra::grid.arrange(grobs = grobs, ncol = ncol)
  } else {
      p
    }
}


#' Plot the variance explained
#'
#' This function first verify that the variance explained has been computed,
#' and if not computes it. It then creates a plot of the variance explained
#' for each PC.
#'
#' @param seur The Seurat object from which the PCA data will be taken
#' @param npcs The number of PCs to plot
#'
#' @return A ggplot object
#' @export
plotVarianceExplained <- function(seur, npcs = 50){
  if (!"variance_explained" %in% names(seur@reductions$pca@misc)){
    seur <- varianceExplained(seur, "pca")
  }
  plot_data <- data.frame(y = seur@reductions$pca@misc$variance_explained[1:npcs], x = 1:npcs)

  ggplot(plot_data, aes(x = x, y = y)) +
    geom_point() +
    labs(title = "Variance explained for each PCs", x = "PC", y = "Variance explained") +
    theme_bw()
}

#' Plot the cumulative variance explained
#'
#' This function first verify that the variance explained has been computed,
#' and if not computes it. It then computes the cumulative sum of the variances,
#' and then creates a plot of the cumulative variance explained for each PC.
#'
#' @param seur The Seurat object from which the PCA data will be taken
#' @param npcs The number of PCs to plot
#'
#' @return A ggplot object
#' @export
plotCumulativeVarExplained <- function(seur, npcs = 50){
  if (!"variance_explained" %in% names(seur@reductions$pca@misc)){
    seur <- varianceExplained(seur, "pca")
  }
  plot_data <- data.frame(y = cumsum(seur@reductions$pca@misc$variance_explained[1:npcs]),
                          x = 1:npcs,
                          variance_explained = seur@reductions$pca@misc$variance_explained)

  ggplot(plot_data, aes(x = x, y = y, color = variance_explained)) +
    geom_point() +
    labs(title = "cumulative variance by PCs", x = "PC", y = "Cumulative variance") +
    theme_bw() +
    ylim(0, NA) +
    scale_color_gradientn(colours = viridis(100, option = "B", end = 0.99, direction = 1), trans = "log10")
}



#' Plot the genes ranked by expression
#'
#' Genes will be ranked by expression, and then plotted, with x axis being the
#' rank and y axis the average expression in the cells, for each ident in group.by.
#' Some genes will then be highlighted.
#'
#' @param seur Seurat object to use.
#' @param genes_highlight List of genes to highlight.
#' If some genes are not found in the dataset, they will just not be plotted.
#' @param group.by Which ident to group the data by. Default to "orig.ident".
#' @param assay Assay to use. Default to "RNA"
#' @param slot slot to take the data from. Using data will assume that the data
#' are log1p normalized, and thus will exp1p them before getting the average expression.
#' Default to "data"
#' @param colmap colmap to use for each ident. Default tu hue_pal()
#' @param nudge_x either a number or "middle". Will nudge the gene name
#' on the right if positive and left if negative.
#' Middle will try to align the gene names in the middle. Default to 2
#' @param ncol number of column to plot, only useful if more than 1 ident.
#' Default to 2
#'
#' @return a ggplot object if only one ident, or a gridExtra object if more than one ident.
#' @export
#'
plotGenesRank <- function(seur, genes_highlight, group.by = "orig.ident", assay = "RNA", slot = "data", colmap = FALSE, nudge_x = "middle", ncol = 2){
  suppressMessages(avg_exp <- as.data.frame(AverageExpression(seur, assays = assay, slot = slot, group.by = group.by)[[assay]]))

  if (is.factor(seur@meta.data[[group.by]])) {
    vars = levels(seur@meta.data[[group.by]])
  } else {
    vars = unique(seur@meta.data[[group.by]])
  }

  if (length(colnames(avg_exp == 1))) {
    colnames(avg_exp) <- vars
  }

  if (colmap == FALSE) {
    colmap = scales::hue_pal()(ncol(avg_exp))
    names(colmap) <- colnames(avg_exp)
  } else if (is.null(names(colmap))) {
    names(colmap) <- colnames(avg_exp)
  }
  avg_exp$names = ifelse(row.names(avg_exp) %in% genes_highlight, row.names(avg_exp), '')

  ps <- lapply(vars, FUN = function(x){
    avg_id <- avg_exp[order(avg_exp[[x]], decreasing = F),]
    avg_id$x <- 1:nrow(avg_id)
    avg_id1 <- subset(avg_id, avg_id$names == '')
    avg_id2 <- subset(avg_id, avg_id$names != '')

    p <- ggplot(avg_id1, aes(x = x, y = .data[[x]])) +
      geom_point(colour = colmap[x], shape = 19) +
      geom_point(data = avg_id2, colour = "black", fill = colmap[x], shape = 21) +
      scale_y_continuous(trans = "log1p") +
      ggtitle(x) +
      labs(y = "Average expression normalized") +
      theme(axis.title.x = element_blank(),
            axis.line = element_line(),
            panel.background=element_blank())

    if (nudge_x == "middle") {
      p <- p + ggrepel::geom_text_repel(data = avg_id2,
                                        aes(label=names),
                                        segment.color = 'grey70',
                                        max.overlaps = 4000,
                                        min.segment.length = 0.1,
                                        point.padding = 0,
                                        force_pull = 0.5,
                                        nudge_x = -max(avg_id$x)%/%2 + max(avg_id$x) - avg_id2$x,
                                        direction = "y")
    } else {
      p <- p + ggrepel::geom_text_repel(data = avg_id2,
                                        aes(label=names),
                                        segment.color = 'grey70',
                                        max.overlaps = 4000,
                                        min.segment.length = 0.1,
                                        point.padding = 0,
                                        force_pull = 0.5,
                                        nudge_x = nudge_x,
                                        direction = "y")
    }
    return(p)
  })
  if (length(ps)>1){
    gridExtra::grid.arrange(grobs = ps, ncol = 2)
  } else {ps[[1]]}
}


