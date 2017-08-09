context("filterClusters")

test_that("valid input parameters", {
  data <- fullDataMapped
  blank <- c()
  brokenData <- data.table::copy(data)
  brokenData[,Collapsed_Peptide := NULL]
  brokenData[,Collapsed_Peptide := NULL]
  cutoff_1 <- 1
  cutoff_0 <- 0
  cutoff_tooBig <- 1.5
  cutoff_tooSmall <- -1
  testthat::expect_error(filterClusters(blank),
                         "The provided input table is not a data.table. Please provide an output table from mapClusters.")
  testthat::expect_error(filterClusters(brokenData),
                         "The provided input table does not contain all required columns. Please provide an output table from mapClusters.")
  testthat::expect_error(filterClusters(data,cutoff_tooBig),
                         "The provided min_main_cluster_freq value is not valid. Please choose a value between 0 and 1.")
  testthat::expect_error(filterClusters(data,cutoff_tooSmall),
                         "The provided min_main_cluster_freq value is not valid. Please choose a value between 0 and 1.")
  testData <- subset(data,Collapsed_Peptide=="SLSGISSSDLTESGALLHDR(2P)")
  testthat::expect_equal(nrow(filterClusters(testData)),196)
  testthat::expect_equal(nrow(filterClusters(testData,cutoff_1)),0)
  testthat::expect_equal(nrow(filterClusters(testData,cutoff_0)),196)
  dataFiltered <- filterClusters(data,0.6)
  testthat::expect_true(all(dataFiltered$new_cluster==dataFiltered$main_cluster_of_tg))
  testthat::expect_gte(min(dataFiltered$main_cluster_tg_freq),0.6)
})
