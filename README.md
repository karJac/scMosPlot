# scMosPlot
simple mosaic plot for single cell data, with option to normalize to cell number per condition
![image](https://github.com/user-attachments/assets/8d8d785f-9289-4b10-8bc2-292221beb0fc)


#' This function creates a mosaic plot based on the provided data frame. <br>
#'<br>
#' @param df_ A data frame with 'cluster_id' and 'condition' columns._<br>
#' @param stat _Optional statistical data, has to be calculated independly using other packages._<br>
#' @param gap_size _Numeric value for the gap size between fills (default: 0.02)._<br>
#' @param spacing _Numeric value for spacing between bars (default: 200)._<br>
#' @param xlab _X-axis label._<br>
#' @param fill_lab _Legend title for the fill._<br>
#' @param title _Plot title._<br>
#' @param alpha_range _Range of alpha transparency values._<br>
#' @param cols _Optional vector of colors._<br>
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
