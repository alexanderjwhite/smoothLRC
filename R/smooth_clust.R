#' Clustering of Right Singular Vectors
#'
#' @inheritParams smooth_lrc
#'
#' @param v matrix; right singular vectors
#'
#' @importFrom mclust Mclust
#'
#' @return vector of cluster labels
#' @export
#'
#' @examples
smooth_clust <- function(v, n_clust){
  Mclust(v, G = n_clust, modelNames = "EEE")$classification
}
