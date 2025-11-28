test_that("is_rec works for simple numeric sequences", {
  Xt <- c(-0.5, -0.4, 0.2, 0.2, -1.0, 0.1, 0.8, -0.2, -0.9, 0.8)
  expect_equal(is_rec(Xt), c(1,1,1,0,0,0,1,0,0,0))
})

test_that("is_rec handles a single-element vector", {
  expect_equal(is_rec(5), 1)
})
