# scMosPlot
simple mosaic plot for single cell data, with option to normalize to cell number per condition

#' This function creates a mosaic plot based on the provided data frame.
#'
#' @param df A data frame with 'cluster_id' and 'condition' columns.
#' @param stat Optional statistical data.
#' @param gap_size Numeric value for the gap size between fills (default: 0.02).
#' @param spacing Numeric value for spacing between bars (default: 200).
#' @param xlab X-axis label.
#' @param fill_lab Legend title for the fill.
#' @param title Plot title.
#' @param alpha_range Range of alpha transparency values.
#' @param cols Optional vector of colors.
#'
#' @return A ggplot object.
#' @export


# How to use:
labels <- c("T cell", "T cell", "T cell", "T cell", "Mph", "Mph","T cell", "T cell", "T cell", "T cell", "Mph", "Mph","T cell", "T cell", "T cell", "T cell", "Mph", "Mph","T cell", "T cell", "T cell", "T cell", "Mph", "Mph")
condition <- c("tumor", "tumor", "tumor", "control", "control", "tumor","tumor", "tumor", "tumor", "control", "control", "tumor","tumor", "tumor", "tumor", "control", "control", "tumor","tumor", "tumor", "tumor", "control", "control", "tumor")
stat <- c(0.1, 0.02)  # statistical significance, one per label type

mosPlot(labels = labels, condition = condition)
mosPlotNorm(labels = labels, condition = condition, stat=stat)


# How to use with Seurat Object v4
mosPlot(labels=seurat_object$labels, condition=seurat_object$condition, stat = c(0.02,0.08,0.2,0.001,0.00002,0.5,0.01)

mosPlotNorm(labels=seurat_object$labels, condition=seurat_object$condition, stat = c(0.02,0.08,0.2,0.001,0.00002,0.5,0.01), cols=rev(cols), plotWidth = T, lineWidth = 500, hjust_lineText = 0.2)
