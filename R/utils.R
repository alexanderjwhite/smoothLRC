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
  logcounts <- as(matrix(log2(counts_vec + 1), nrow = n_gene), "dgCMatrix")

  col_data <- list()
  spatial_grid <- expand.grid(row = seq(1, n_row), col = seq(1, n_col))
  col_data$row <- spatial_grid$row
  col_data$col <- spatial_grid$col
  col_data <- as.data.frame(do.call(cbind, col_data))
  sce <- SummarizedExperiment(assays=list(counts=counts,logcounts=logcounts), colData=col_data)
  rownames(sce) <- paste0("Feature ", seq(1, nrow(sce)))
  colnames(sce) <- paste0("Pixel ", seq(1, ncol(sce)))
  sce
}
