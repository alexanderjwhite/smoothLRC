#' smoothLRC algorithm
#'
#' @inheritParams smooth_lrc
#' @param u_init matrix; u initialization matrix
#' @param v_init matrix; v initialization matrix
#' @param w matrix; distance matrix. If null, computed to specification with k nearest neighbors.
#' @param index w indicies.
#'
#' @return SummarizedExperiment object with u, v and cluster labels.
#' @export
#'
#' @examples
smooth_model <- function(x, u_init, v_init, w, index, lambda, k, epsilon, maxiter){

  model <- fct_c_optimize(x = input, u = u_init, v = v_init, w = w, index = index, lambda = lambda, epsilon = epsilon, maxiter = maxiter)

  return(list(u = model$u, v = model$v))

}
