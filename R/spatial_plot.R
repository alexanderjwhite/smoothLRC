#' Create Spatial Plot of Feature Expression or Cluster Labels
#'
#' @inheritParams smooth_lrc
#' @param type Plot type, "raw", "smooth, or "cluster". If "raw", the expression of feature_name using assay_type is plotted. If "smooth" the smoothed expression of feature_name is plotted. If "cluster" the cluster labels are plotted. If type="cluster" and labels=NULL, the smoothLRC labels are used.
#' @param feature_name If type="raw" or type="smooth", the name of the feature to be plotted.
#' @param labels If type="cluster" a vector of labels for each observation. This should be the same length as the number of columns of input. If NULL, smoothLRC cluster labels are used.
#' @param assay_type If type="raw" or type="smooth", the SummarizedExperiment assay type.
#' @param spot_size The size of each pixel.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' sce <- example_sce()
#' spatial_plot(sce, feature_name = "Feature 1", spot_size = 10)
spatial_plot <- function(input, type = c("raw", "smooth", "cluster"), feature_name = NULL, labels = NULL, assay_type = "logcounts", spot_size = 1.5){

  if(length(type) > 1){
    type <- type[1]
  }

  if(type != "cluster" & is.null(feature_name)){
    stop("A feature must be provided for an expression plot.")
  }

  title <- ""
  feature_row <- which(rownames(input) == feature_name)
  plot_data <- as.data.frame(SummarizedExperiment::colData(input)[,c("row", "col")])

  if(type != "cluster"){

    if(type == "smooth"){
      u <- input@metadata$smooth_u
      v <- t(input@metadata$smooth_v)
      smoothed <- u%*%v
      feature_data <- smoothed[feature_row,]
      title <- paste("Smoothed", feature_name)
    } else {
      feature_data <- SummarizedExperiment::assay(input, assay_type)[feature_name, ]
      title <- paste("Raw", assay_type, feature_name)
    }

    plot_data <- cbind(plot_data, data.frame(feature = feature_data))
    plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = col, y = -row, color = feature)) +
      ggplot2::geom_point(size = spot_size) +
      ggplot2::scale_color_viridis_c() +
      ggplot2::ggtitle(title)

  } else {
    smoothlrc_labels <- input@colData$smooth_cluster
    if(is.null(smoothlrc_labels) & is.null(labels)){
      stop("For cluster plot, labels are not provided and smoothLRC cluster labels do not exist.")
    } else if(is.null(labels)){
      labels <- smoothlrc_labels
      title <- "smoothLRC Clusters"
    }
    plot_data <- cbind(plot_data, data.frame(cluster = labels))
    plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = col, y = -row, color = cluster)) +
      ggplot2::geom_point(size = spot_size) +
      ggplot2::scale_color_viridis_d() +
      ggplot2::ggtitle(title)
  }

  plot +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.title = ggplot2::element_blank(),
                   axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   legend.title = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank())



}
