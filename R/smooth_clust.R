#' Clustering of Right Singular Vectors
#'
#' @inheritParams smooth_lrc
#'
#' @param v matrix; right singular vectors
#'
#' @importFrom mclust Mclust mclustBIC
#'
#' @return vector of cluster labels
#' @export
#'
smooth_clust <- function(v, n_clust){
  Mclust(v, G = n_clust, modelNames = "EEE")$classification
}
