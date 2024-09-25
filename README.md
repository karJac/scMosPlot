# scMosPlot
simple mosaic plot for single cell data, with option to normalize to cell number per condition
![image](https://github.com/user-attachments/assets/8d8d785f-9289-4b10-8bc2-292221beb0fc)


#' This function creates a mosaic plot based on the provided data frame. <br>
#'<br>
#' @param df A data frame with 'cluster_id' and 'condition' columns.<br>
#' @param stat Optional statistical data.<br>
#' @param gap_size Numeric value for the gap size between fills (default: 0.02).<br>
#' @param spacing Numeric value for spacing between bars (default: 200).<br>
#' @param xlab X-axis label.<br>
#' @param fill_lab Legend title for the fill.<br>
#' @param title Plot title.<br>
#' @param alpha_range Range of alpha transparency values.<br>
#' @param cols Optional vector of colors.<br>
#'<br>
#' @return A ggplot object.<br>
#' @export<br>


# How to use:
labels <- c("T cell", "T cell", "T cell", "T cell", "Mph", "Mph","T cell", "T cell", "T cell", "T cell", "Mph", "Mph","T cell", "T cell", "T cell", "T cell", "Mph", "Mph","T cell", "T cell", "T cell", "T cell", "Mph", "Mph")<br>
condition <- c("tumor", "tumor", "tumor", "control", "control", "tumor","tumor", "tumor", "tumor", "control", "control", "tumor","tumor", "tumor", "tumor", "control", "control", "tumor","tumor", "tumor", "tumor", "control", "control", "tumor")<br>
stat <- c(0.1, 0.02)  # statistical significance, one per label type<br>

mosPlot(labels = labels, condition = condition)<br>
mosPlotNorm(labels = labels, condition = condition, stat=stat)<br>


# How to use with Seurat Object v4
mosPlot(labels=seurat_object$labels, condition=seurat_object$condition, stat = c(0.02,0.08,0.2,0.001,0.00002,0.5,0.01)<br>
<br>
mosPlotNorm(labels=seurat_object$labels, condition=seurat_object$condition, stat = c(0.02,0.08,0.2,0.001,0.00002,0.5,0.01), cols=rev(cols), plotWidth = T, lineWidth = 500, hjust_lineText = 0.2)<br>
