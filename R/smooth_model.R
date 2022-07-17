#' smoothLRC algorithm
#'
#' @inheritParams smooth_lrc
#' @param x sparse matrix; assay matrix
#' @param u_init matrix; u initialization matrix
#' @param v_init matrix; v initialization matrix
#' @param w matrix; distance matrix. If null, computed to specification with k nearest neighbors.
#' @param index w indicies.
#'
#' @return SummarizedExperiment object with u, v and cluster labels.
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
#' smooth_model(x, init$u0, init$v0, init$w, init$index, lambda, epsilon, maxiter)
#'
smooth_model <- function(x, u_init, v_init, w, index, lambda, epsilon, maxiter){

  model <- fct_c_optimize(x, u_init, v_init, w, index, lambda, epsilon, maxiter)

  return(list(u = model$u, v = model$v))

}
