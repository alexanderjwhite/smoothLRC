## ========================
## Setup
## ========================
sce <- example_sce()
x <- SummarizedExperiment::assay(sce)
coords <- SummarizedExperiment::colData(sce)[, c("col", "row")]
k <- 10
knn <- 6
init <- smooth_init(x, coords, k, knn)
lambda <- 1
epsilon <- 1e-3
maxiter <- 5

## ========================
## Test smooth_init
## ========================

test_that("model returns list with appropriate attributes", {
  model <- smooth_model(x, init$u0, init$v0, init$w, init$index, lambda, epsilon, maxiter)

  expect_type(model, "list")
  expect_true("u" %in% names(model))
  expect_true("v" %in% names(model))
  expect_true(nrow(model$u) == nrow(x))
  expect_true(ncol(model$u) == k)
  expect_true(nrow(model$v) == ncol(x))
  expect_true(ncol(model$v) == k)

})
