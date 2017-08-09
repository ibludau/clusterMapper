context("generateData")

test_that("valid input parameters", {
  data <- exampleData
  blank <- c()
  brokenData <- data.table::copy(data)
  brokenData[,Collapsed_Peptide := NULL]
  out <- filterClusters(mapClusters(filterData(exampleData)))
  out_false <- exampleData
  out_broken <- data.table::copy(out)
  out_broken[, new_cluster := NULL]
  finalOut <- generateMatrix(data,out)
  testthat::expect_error(generateMatrix(blank,out),
                         "The provided input table \"rawInput\" is not a data.table.")
  testthat::expect_error(generateMatrix(brokenData,out),
                         "The provided input table \"rawInput\" does not contain all required columns.")
  testthat::expect_error(generateMatrix(data,blank),
                         "The provided input table \"processedData\" is not a data.table.")
  testthat::expect_error(generateMatrix(data,out_false),
                         "The provided input table \"processedData\" does not contain all required columns.")
  testthat::expect_error(generateMatrix(data),
                         "argument \"processedData\" is missing, with no default")
  testthat::expect_error(generateMatrix(processedData = data),
                         "argument \"rawInput\" is missing, with no default")
  testthat::expect_equal(nrow(out), nrow(finalOut))
  testthat::expect_gte(nrow(data), nrow(finalOut))
  testthat::expect_true(all(c("transition_group_id","run_id","RT","Collapsed_Peptide","Cluster","new_cluster") %in% names(finalOut)))
 })
