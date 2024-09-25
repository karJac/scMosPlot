#' mymosplot Function
#'
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

# Load necessary libraries
library(ggplot2)
library(ggmosaic)
library(scales)

# create mosplot graph function
mymosplot <- function(df, stat = NULL, gap_size = 0.02, spacing = 200, xlab = "Population", 
                      fill_lab = "Condition", title = NULL, alpha_range = c(0.3, 1), cols = NULL) {
  if (is.null(cols)) {
    cols <- rainbow(length(levels(df$condition)), s = 0.8, v = 0.8)
  }

  # Calculate total counts per cluster_id
  cluster_counts <- df %>%
    group_by(cluster_id) %>%
    summarise(total_count = n())
  if (!is.null(stat)) {
    alpha_levels <- ifelse(stat > 0.05, 0.4, 1)
    cluster_counts$alpha <- alpha_levels[match(cluster_counts$cluster_id, names(alpha_levels))]
  } else {
    cluster_counts$alpha <- 1
  }
  cluster_counts$alpha <- as.factor(cluster_counts$alpha)

  # Calculate counts per cluster_id and condition
  cluster_condition_counts <- df %>%
    group_by(cluster_id, condition) %>%
    summarise(count = n()) %>%
    ungroup()

  # Set spacing between bars
  spacing <- spacing # Adjust this value to change the distance between columns
  gap_size <- gap_size # Adjust this value to change the gap between fills (0-1)

  # Compute bar positions
  cluster_counts <- cluster_counts %>%
    mutate(width = total_count)

  # Initialize positions
  cluster_counts$xmin <- NA
  cluster_counts$xmax <- NA

  cluster_counts$xmin[1] <- 0
  cluster_counts$xmax[1] <- cluster_counts$width[1]

  if (nrow(cluster_counts) >= 2) {
    for (i in 2:nrow(cluster_counts)) {
      cluster_counts$xmin[i] <- cluster_counts$xmax[i - 1] + spacing
      cluster_counts$xmax[i] <- cluster_counts$xmin[i] + cluster_counts$width[i]
    }
  }

  # Merge data for plotting
  data_for_plot <- cluster_condition_counts %>%
    left_join(cluster_counts %>% dplyr::select(cluster_id, xmin, xmax, total_count, alpha), by = "cluster_id")

  # Compute stacking positions
  data_for_plot2 <- data_for_plot %>%
    group_by(cluster_id) %>%
    arrange(condition) %>% # Adjust if you have a specific order for conditions
    mutate(
      prop = count / sum(count), # Calculate the proportion
      n_segments = n(),
      total_gap = gap_size * (n_segments - 1),
      adjusted_prop = prop * (1 - total_gap),
      cumulative_adjusted_prop = cumsum(adjusted_prop),
      cumulative_gap = gap_size * (row_number() - 1),
      ymax = cumulative_adjusted_prop + cumulative_gap,
      ymin = ymax - adjusted_prop
    )

  
  # Create the plot
  p_my <- ggplot(data_for_plot2) +
    geom_rect(
      aes(
        xmin = xmin, xmax = xmax,
        ymin = ymin, ymax = ymax,
        fill = condition, alpha = alpha
      ),
      color = "black"
    ) +
    scale_x_continuous(
      breaks = (cluster_counts$xmin + cluster_counts$xmax) / 2,
      labels = cluster_counts$cluster_id,
      expand = c(0, 0)
    ) +
    scale_y_continuous(expand = c(0, 0),
                       breaks = c(0.0, 0.25, 0.5, 0.75, 1.0)) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 14), # X-axis labels
      axis.text.y = element_text(size = 11), # Y-axis labels
      axis.title.x = element_text(size = 16), # X-axis title
      axis.title.y = element_text(size = 16), # Y-axis title
      legend.text = element_text(size = 13), # Legend text
      legend.title = element_text(size = 15), # Legend title
      plot.title = element_text(size = 16), # Plot title
      axis.line.y = element_blank() # Remove Y axis line
    ) +
    labs(x = xlab, y = "Proportion", title = title) +
    scale_fill_manual(
      values = cols,
      breaks = rev(levels(data_for_plot2$condition)),
      name = fill_lab # Use your existing legend title variable
    ) +
    scale_alpha_discrete(range = alpha_range, labels = c("> 0.05", "<= 0.05")) +
    guides(alpha = guide_legend(title = "pval_adj"), fill = guide_legend(title = fill_lab))

  if (length(unique(data_for_plot2$alpha)) == 1) { # otherwise the code was giving minimal alpha if you dont provide stat (then $alpha == 1)
    if (unique(data_for_plot2$alpha) == 1) {
      p_my <- p_my +
        scale_alpha_discrete(range = c(0.99, 1)) +
        guides(alpha = FALSE)
    }
  }

  return(p_my)
}




############################################






# mosPlot function
mosPlot <- function(labels, condition, stat = NULL, cols = NULL, legend = TRUE, 
                    gap_size=0.02, spacing=200, alpha_range = c(0.3,1),
                    title = 'Percentage of Cells from Different Conditions',
                    x_lab = "Population", fill_lab = "Condition", 
                    middleLine = FALSE, plotWidth = TRUE, lineWidth = 1000, hjust_lineText = 0.2) {
  
  # Create data frame
  x <- data.frame(cluster_id = labels, condition = condition)
  
  # Convert to factors with levels based on unique values if is not already a factor
  if (!is.factor(x$cluster_id)){
  x$cluster_id <- factor(x$cluster_id, levels = unique(labels))
  }
  
  if (!is.factor(x$condition)){
  x$condition <- factor(x$condition, levels = unique(condition))
  }
  
  # Set default colors if not provided
  if (is.null(cols)) {
    cols <- rainbow(length(levels(x$condition)),s=0.8,v=0.8)
  }
  
  # If stat is provided, calculate alpha based on statistical significance
  if (!is.null(stat)) {
    # Ensure stat is a named vector matching cluster IDs
    if (is.null(names(stat))) {
      unique_clusters <- levels(x$cluster_id)
      if (length(stat) != length(unique_clusters)) {
        stop("Length of stat must be equal to the number of clusters")
      }
      names(stat) <- unique_clusters
    }
    
    # Map stat values to alpha levels
    alpha_levels <- ifelse(stat > 0.05, 0.4, 1)
    
    # Add alpha column to x
    x$alpha <- alpha_levels[as.character(x$cluster_id)] 
  } else {
    # If stat is not provided, set alpha to 1
    x$alpha <- 1
  }
  
  # Base plotp 
  p <- mymosplot(df = x, stat = stat, cols = cols, gap_size = gap_size, spacing = spacing, fill_lab = fill_lab, xlab = xlab, title = title)
  
  if (legend == FALSE){
    p <- p + NoLegend()
  }
    
  # Add additional elements if 'middleLine' is TRUE
  if (middleLine) {
    p <- p +
      geom_hline(yintercept = 0.5)
  }
  
  #add line showing width of 1000cells
  if (plotWidth) {
    p <- p +
      geom_segment(
        aes(x = 0, xend = lineWidth, y = 1.05, yend = 1.05), # Adjust y position above the plot
        color = "black", size = 1.2
      ) +
      annotate(
        "text", x = lineWidth / 2, y = 1.08, label = paste0(as.character(lineWidth), " cells"), size = 5, hjust = hjust_lineText, vjust=0.3
      ) + 
      theme(plot.title = element_text(vjust = 1.5)) +
      coord_cartesian(ylim = c(0, 1.15))
  }
  
  return(p)
}



# mosPlotNorm function
mosPlotNorm <- function(labels, condition, stat = NULL, cols = NULL, legend = TRUE, 
                        gap_size=0.02, spacing=200, alpha_range = c(0.3,1),
                        title = 'Percentage of Cells from Different Conditions, Normalized by Condition',
                        x_lab = "Population", fill_lab = "Condition", 
                        middleLine = FALSE, plotWidth = TRUE, lineWidth = 1000, hjust_lineText = 0.2) {
  
  # Create data frame
  x <- data.frame(cluster_id = labels, condition = condition)
  
  # Create contingency table
  tab <- table(x$cluster_id, x$condition)
  
  # Calculate total counts per condition
  total_counts <- colSums(tab)
  
  # Calculate scaling factors to equalize total counts per condition
  scaling_factors <- mean(total_counts) / total_counts
  
  # Apply scaling factors
  norm_tab <- sweep(tab, 2, scaling_factors, FUN = "*")
  
  # Round normalized counts to integers
  norm_tab <- round(norm_tab)
  
  # Reconstruct data frame with normalized counts
  norm_x <- as.data.frame(as.table(norm_tab))
  colnames(norm_x) <- c("cluster_id", "condition", "Freq")
  
  # Repeat rows according to 'Freq'
  norm_x <- norm_x[rep(seq_len(nrow(norm_x)), norm_x$Freq), 1:2]
  
  
  # Call mosPlot with normalized data
  p <- mosPlot(labels = norm_x$cluster_id, condition = norm_x$condition, legend = TRUE,
              cols = cols, stat = stat, middleLine = middleLine,
               title = title, plotWidth = plotWidth, lineWidth = lineWidth, hjust_lineText = hjust_lineText)
  
  return(p)
}
