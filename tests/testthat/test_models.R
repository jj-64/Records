context("basic model behavior")

test_that("IID mean/var are consistent", {
  expect_equal(rec_count_mean_iid(5), sum(1/(1:5)))
  expect_true(rec_count_var_iid(5) >= 0)
})

test_that("DTRW mean approx works", {
  expect_true(rec_count_mean_DTRW(50, approx = TRUE) > 0)
})

test_that("dist returns probability-like values", {
  p <- rec_count_dist_DTRW(m = 3, T = 10, approx = TRUE)
  expect_true(is.numeric(p))
  expect_true(p >= 0)
})