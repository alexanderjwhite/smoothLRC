#' Initialization for smoothLRC
#'
#' @inheritParams smooth_lrc
#' @param x matrix; count assay.
#' @param coords matrix; spatial coordinates.
#' @param knn integer; nearest neighbors.
#'
#' @importFrom sparsesvd sparsesvd
#' @importFrom methods as
#' @importFrom FNN get.knn
#' @importFrom Matrix sparseMatrix
#'
#' @return list of objects for smooth_model function.
#' @export
#'
#' @examples
#'
#' sce <- example_sce()
#' x <- SummarizedExperiment::assay(sce)
#' coords <- SummarizedExperiment::colData(sce)[, c("col", "row")]
#' k <- 10
#' smooth_init(x, coords, k)
#'
smooth_init <- function(x, coords, k, knn = 6){

  sp_svd <- sparsesvd::sparsesvd(log(x+1),rank = k)
  u0 <- sp_svd$u %*% sqrt(diag(sp_svd$d))
  v0 <- sp_svd$v %*% sqrt(diag(sp_svd$d))

  nn <- FNN::get.knn(coords, k = knn)
  index <- nn$nn.index
  neighbours <- as.vector(t(index))
  w <- Matrix::sparseMatrix(i = rep(seq(1, nrow(coords)), each = knn), j = neighbours, x = 1, dims = c(nrow(coords), nrow(coords)))
  w <- methods::as(w, "dgCMatrix")

  return(list(u0 = u0, v0 = v0, w = w, index = index))

}
