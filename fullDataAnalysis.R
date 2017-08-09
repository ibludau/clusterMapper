library(devtools)
library(data.table)

setwd("/IMSB/ra//ibludau/clusterMapper/clusterMapper/")
load_all()

data <- fread("/IMSB/ra/ibludau/phospho_cluster_mapping/input/combined_clusters_head_all_mscore.txt")

dataMapped <- fullDataMapped
dataMappedFiltered <- filterClusters(dataMapped)
dataFinal <- generateMatrix(rawInput = data, processedData = dataMappedFiltered)
