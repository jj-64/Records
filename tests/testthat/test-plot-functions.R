```r
context("Plot functions")

test_that("Lorenz curve runs without error", {
  data(SumData)
  expect_silent(plot_lorenz(SumData[1, ]))
})

test_that("Decile plot runs without error", {
  data(SumData)
  expect_silent(plot_decile_shares(SumData[1, ]))
})

test_that("Headcount model plot runs without error", {
  data(HC_limited_data)
  data(CI_limited_data)
  expect_silent(plot_headcount_models(HC_limited_data[1, ], CI_limited_data[CI_limited_data$Country==HC_limited_data$Country[1], ]))
})

test_that("Country summary plot works", {
  data(SumData)
  skip_if_not("patchwork" %in% rownames(installed.packages()))
  expect_silent(plot_country_summary(SumData$Country[1]))
})
```
