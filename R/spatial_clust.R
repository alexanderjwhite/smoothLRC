#' Main Spatial Clustering Method
#'
#' @param input SingleCellExperiment object containing counts assay and row/col coordinates.
#' @param u_init matrix; u initialization matrix
#' @param v_init matrix; v initialization matrix
#' @param w matrix; distance matrix. If null, computed to specification with k nearest neighbors. 
#' @param index w index.
#' @param lambda numeric; penalization parameter. set to NULL for lambda selection.
#' @param cnorm a numeric value for the amount of V normalization. NULL for no normalization.
#' @param k integer; if knn = TRUE, number of nearest neighbors to consider
#' @param epsilon numeric; convergence criterion
#' @param maxiter integer; maximum number of iterations
#' 
#' @importFrom methods as
#' @importFrom SummarizedExperiment assay colData
#' @importFrom FNN get.knn
#' @importFrom Matrix sparseMatrix

#'
#' @return Object containing u, v and convergence information
#' @export
#'
#' @examples 
spatial_clust <- function(input, u_init, v_init, w = NULL, index = NULL, lambda = 1, cnorm = NULL, k = 20, epsilon = 1e-3, maxiter = 1e3){
  
  # x <- SummarizedExperiment::assay(input)
  
  if(is.null(w)){
    coords <- SummarizedExperiment::colData(input)[, c("row", "col")]
    nn <- FNN::get.knn(coords, k = k)
    index <- nn$nn.index
    neighbours <- as.vector(t(index))
    w <- Matrix::sparseMatrix(i = rep(1:nrow(coords), each = k), j = neighbours, x = 1, dims = c(nrow(coords), nrow(coords)))
    w <- methods::as(w, "dgCMatrix")
  }
  
  fct_c_optimize(x = input, u = u_init, v = v_init, w = w, index = index, lambda = lambda, cnorm = cnorm, epsilon = epsilon, maxiter = maxiter)

}