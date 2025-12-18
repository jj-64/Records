
test_that("rec_count returns correct number of records", {
  Xt <- c(-0.5, -0.4, 0.2, 0.2, -1.0, 0.1, 0.8, -0.2, -0.9, 0.8)
  expect_equal(rec_count(Xt), 4)
})
