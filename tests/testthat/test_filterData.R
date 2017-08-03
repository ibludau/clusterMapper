context("filterData")

test_that("valid input parameters", {
  data <- exampleData
  blank <- c()
  brokenData <- data.table::copy(data)
  brokenData[,Collapsed_Peptide := NULL]
  cutoff_1 <- 1
  cutoff_0 <- 0
  cutoff_tooBig <- 1.5
  cutoff_tooSmall <- -1
  testthat::expect_error(filterData(blank),
                         "The provided input table is not a data.table.")
  testthat::expect_error(filterData(brokenData),
                         "The provided input table does not contain all required columns.")
  testthat::expect_error(filterData(data,cutoff_tooBig),
                         "The provided min_run_coverage value is not valid. Please choose a value between 0 and 1.")
  testthat::expect_error(filterData(data,cutoff_tooSmall),
                         "The provided min_run_coverage value is not valid. Please choose a value between 0 and 1.")
  testthat::expect_equal(nrow(filterData(data)),56)
  testthat::expect_equal(nrow(filterData(data,cutoff_1)),42)
  testthat::expect_equal(nrow(filterData(data,cutoff_0)),nrow(data))
})
