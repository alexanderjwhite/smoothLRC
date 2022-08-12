#' Cluster the right singular vectors from model output
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
#' @examples
#'
#' sce <- example_sce()
#' x <- SummarizedExperiment::assay(sce)
#' coords <- SummarizedExperiment::colData(sce)[, c("col", "row")]
#' k <- 10
#' init <- smooth_init(x, coords, k)
#' lambda <- 1
#' epsilon <- 1e-3
#' maxiter <- 5
#' model <- smooth_model(x, init$u0, init$v0, init$w, init$index, lambda, epsilon, maxiter)
#' n_clust <- 3
#' smooth_clust(model$v, n_clust)
#'
smooth_clust <- function(v, n_clust){
  Mclust(v, G = n_clust, modelNames = "EEE")$classification
}
