## ========================
## Test utils
## ========================

test_that("example summarized experiment object is created", {
  sce <- example_sce()
  expect_type(sce, "S4")
})
