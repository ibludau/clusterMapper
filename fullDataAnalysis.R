library(devtools)
library(data.table)

#setwd("/IMSB/ra//ibludau/clusterMapper/clusterMapper/")
setwd("/Volumes/ibludau-1/clusterMapper/clusterMapper/")
load_all()

#data <- fread("/IMSB/ra/ibludau/phospho_cluster_mapping/input/combined_clusters_head_all_mscore.txt")
data <- fread("/Volumes/ibludau-1/phospho_cluster_mapping/input/combined_clusters_head_all_mscore.txt")
dataFiltered <- filterData(data,
                           min_run_coverage = 0.6)
dataMapped <- mapClusters(dataFiltered)
#dataMapped <- fullDataMapped
dataMappedFiltered <- filterClusters(dataMapped,
                                     min_main_cluster_freq = 0.6)
dataFinal <- generateMatrix(rawInput = data, processedData = dataMappedFiltered)

write.table(dataFinal,"/Volumes/ibludau-1/clusterMapper/phosphoQTL/dataFinal.tsv", quote = F, sep="\t", row.names = F)
save.image("/Volumes/ibludau-1/clusterMapper/phosphoQTL/workspace_data.RData")
