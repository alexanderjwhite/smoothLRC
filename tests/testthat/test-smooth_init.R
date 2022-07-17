## ========================
## Setup
## ========================
sce <- example_sce()
x <- SummarizedExperiment::assay(sce)
coords <- SummarizedExperiment::colData(sce)[, c("col", "row")]
k <- 10
knn <- 6

## ========================
## Test smooth_init
## ========================

test_that("initialization returns list with appropriate attributes", {
  init <- smooth_init(x, coords, k, knn)

  expect_type(init, "list")
  expect_true("u0" %in% names(init))
  expect_true("v0" %in% names(init))
  expect_true("w" %in% names(init))
  expect_true("index" %in% names(init))
  expect_true(nrow(init$u0) == nrow(x))
  expect_true(ncol(init$u0) == k)
  expect_true(nrow(init$v0) == ncol(x))
  expect_true(ncol(init$v0) == k)
  expect_true(nrow(init$w) == ncol(x))
  expect_true(ncol(init$w) == ncol(x))
  expect_true(nrow(init$index) == ncol(x))
  expect_true(ncol(init$index) == knn)
})
