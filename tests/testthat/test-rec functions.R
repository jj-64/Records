test_that("rec functions handle trivial input correctly", {
  expect_equal(rec_counts(10), 1)
  expect_equal(rec_values(10), 10)
  expect_equal(rec_times(10), 1)
  expect_equal(rec_gaps(10), numeric(0))
})
