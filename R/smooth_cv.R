#' Parameter Selection via Cross-Validation
#'
#' @inheritParams smooth_lrc
#' @param test_size The number of pixels to use in the test set.
#'
#' @return A list containing the penalized likelihood from the model and the cross-validated likelihood from the test set.
#' @export
#'
#' @examples
#' sce <- example_sce()
#' smooth_cv(sce, 1, 10, 100, maxiter = 10)
smooth_cv <- function(input, lambda, k, test_size, epsilon = 1e-3, maxiter = 1e3){

  # Get assay/coordinates from SummarizedExperiment object
  x <- assay(input)
  coords <- colData(input)[, c("col", "row")]

  # Initialize with full input
  print("Initializing components...")
  full_nn <- FNN::get.knn(coords, k = 6)
  full_index <- full_nn$nn.index
  full_dist <- full_nn$nn.dist

  # Find non-neighboring pixels
  test_cols <- NULL
  neighbor_map <- NULL
  n <- 0
  iter <- 0
  while ((n <= test_size) & iter < nrow(coords)){
    sample_i <- sample(1:nrow(coords), 1, replace = FALSE)
    if ((!(sample_i %in% neighbor_map)) & (sum(full_index[sample_i,] %in% test_cols) == 0)){
      test_cols <- c(test_cols,sample_i)
      neighbor_map <- c(neighbor_map,sample_i,full_index[sample_i,])
      n <- n + 1
    }
    iter <- iter + 1
  }

  # Separate into train and test sets
  test_x <- x[,test_cols]
  test_index <- full_index[test_cols,]
  index_map <- matrix(matrix(1:nrow(coords), ncol = 1)[-test_cols,],ncol= 1)
  train_x <- x[,-test_cols]
  train_coords <- coords[-test_cols,]

  # Initialize training components
  init_model <- smooth_init(train_x, train_coords, k)

  # Run model
  print("Running smoothLRC on training set...")
  model <- smooth_model(train_x, init_model$u0, init_model$v0, init_model$w, init_model$index, lambda, epsilon, maxiter)

  # Compute neighborhood likelihood
  print("Compute neighborhood likelihood...")
  cv_like <- mean(as.vector(obs_log_like(test_x, model$u, model$v, test_index, index_map)))
  print("Done!")

  return(list(penal_like = model$penal_lik, cv_like = cv_like))
}
