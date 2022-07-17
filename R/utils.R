#' Load Example sce dataset
#'
#' @param n_row integer; number of simulated spatial pixel rows
#' @param n_col integer; number of simulated spatial pixel columns
#' @param n_gene integer; number of simulated genes
#' @param lambda numeric; poisson parameter
#' @param frac_sparse numeric; fraction of 0s.
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom stats rpois
#' @importFrom methods as
#'
#' @return SummarizedExperiment object
#' @export
#'
#' @examples
#'
#' example_sce()
#'
example_sce <- function(n_row = 10, n_col = 10, n_gene = 1000, lambda = 2, frac_sparse = 0.95){

  n_pixel <- n_row*n_col
  nonzero <- (1-frac_sparse)*n_gene*n_pixel
  nonzero_counts <- rpois(nonzero, lambda)
  counts_vec <- sample(c(rep(0,frac_sparse*n_pixel*n_gene),nonzero_counts), n_pixel*n_gene)
  counts <- as(matrix(counts_vec, nrow = n_gene), "dgCMatrix")

  col_data <- list()
  col_data$row <- rep(seq(1, n_row), each = n_col)
  col_data$col <- rep(seq(1, n_col), each = n_row)
  col_data <- as.data.frame(do.call(cbind, col_data))
  SummarizedExperiment(assays=list(counts=counts), colData=col_data)

}
