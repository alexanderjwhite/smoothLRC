#' Compute the Cross Validated Log Likelihood
#'
#' @param test_x test matrix
#' @param u model computed u
#' @param v model computed v
#' @param test_nn nearest neighbor matrix
#' @param index_map index map
#'
#' @return double
#' @export
#'
#' @examples
cv_like <- function(test_x, u, v, test_nn, index_map){
  obs_log_like(test_x, u, v, test_nn, index_map)
}