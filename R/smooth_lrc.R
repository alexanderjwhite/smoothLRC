#' Cluster spatial transcriptomics data
#'
#' @param input SummarizedExperiment object containing counts assay and row/col coordinates.
#' @param lambda positive numeric; penalization parameter.
#' @param k integer; indicates desired rank of singular value decomposition.
#' @param n_clust integer; number of clusters.
#' @param epsilon positive numeric; convergence criterion.
#' @param maxiter positive integer; maximum desired iterations
#'
#' @importFrom SummarizedExperiment assay colData
#'
#' @return SummarizedExperiment object with u, v and cluster labels.
#' @export
#'
#' @examples
#'
#' sce <- example_sce()
#' lambda <- 5
#' k <- 10
#' n_clust <- 5
#' epsilon <- 1e-3
#' maxiter <- 5
#' smooth_lrc(sce, lambda, k, n_clust, epsilon, maxiter)
#'
smooth_lrc <- function(input, lambda, k, n_clust, epsilon = 1e-3, maxiter = 1e3){

  # Get assay/coordinates from SummarizedExperiment object
  x <- assay(input)
  coords <- colData(input)[, c("col", "row")]

  # Initialize components
  print("Initializing components...")
  init_model <- smooth_init(x, coords, k)

  # Run the model
  print("Running smoothLRC...")
  model <- smooth_model(x, init_model$u0, init_model$v0, init_model$w, init_model$index, lambda, epsilon, maxiter)

  # Cluster the right singular vectors
  print("Clustering right singular vectors...")
  clust <- smooth_clust(model$v, n_clust)

  # Attach components and return
  input@colData$smooth_cluster <- as.factor(unname(clust))

  input@metadata$smooth_u <- model$u
  input@metadata$smooth_v <- model$v
  print("Done!")

  return(input)

}
