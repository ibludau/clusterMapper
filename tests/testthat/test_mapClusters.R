context("mapClusters")

test_that("one cluster", {
  # test easy dataset
  testData <- subset(exampleDataFiltered,Collapsed_Peptide=="KPSLEFLPQPASSTNLNFNSNK(1P)")
  testthat::expect_equal(unique(mapClusters(testData)$new_cluster),1)
  rm(testData)
  # test data with one outlier
  testData <- subset(exampleDataFiltered,Collapsed_Peptide=="NLLTNTPVDYEYESDAEEEQGDKDIK(1P)")
  testMap <- mapClusters(testData)
  outlier_idx <- which((testMap$transition_group_id=="9994_NLLTNTPVDY(UniMod:21)EYESDAEEEQGDKDIK_3_run0") & (testMap$run_id=="0_1"))
  testthat::expect_equal(testMap[outlier_idx]$new_cluster,0)
  testthat::expect_equal(unique(testMap[-outlier_idx]$new_cluster),1)
  rm(testData,testMap,outlier_idx)
  # test data with multiple clusters
  testData <- complexExample
  testMap <- mapClusters(testData)
  testthat::expect_equal(length(which(testMap$new_cluster == 0)),7)
  testthat::expect_equal(length(which(testMap$new_cluster == 1)),85)
  testthat::expect_equal(length(which(testMap$new_cluster == 2)),0)
  testthat::expect_equal(length(which(testMap$new_cluster == 3)),26)
  rm(testData,testMap)
  # test data that previously made problems
  testData <- subset(fullDataFiltered,Collapsed_Peptide == "NPTKSPPPPPSPSTMDTGTSNSPSK(3P)")
  testMap <- mapClusters(testData)
  testthat::expect_equal(length(which(testMap$new_cluster == 0)),1049)
  testthat::expect_equal(length(which(testMap$new_cluster == 1)),0)
  testthat::expect_equal(length(which(testMap$new_cluster == 2)),116)
  testthat::expect_equal(length(which(testMap$new_cluster == 3)),1480)
  rm(testData,testMap)
  # and another problematic one
  testData <- subset(fullDataFiltered,Collapsed_Peptide == "NSNNSFLNSVPHSVTR(2P)")
  testMap <- mapClusters(testData)
  testthat::expect_equal(length(which(testMap$new_cluster == 0)),115)
  testthat::expect_equal(length(which(testMap$new_cluster == 1)),275)
  rm(testData,testMap)
})
