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
smooth_model <- function(x, u_init, v_init, w, index, lambda, epsilon, maxiter){

  model <- fct_c_optimize(x, u_init, v_init, w, index, lambda, epsilon, maxiter)

  return(list(u = model$u, v = model$v))

}
