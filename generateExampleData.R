library('data.table')
data <- fread("/IMSB/ra/ibludau/phospho_cluster_mapping/input/combined_clusters_head_all_mscore.txt")
dataSub <- subset(data,run_id %in% c("0_0","0_1","0_2","0_3","0_4","0_5"))
dataSub <- subset(dataSub,Collapsed_Peptide %in% unique(dataSub$Collapsed_Peptide)[1:5])
exampleData <- copy(dataSub)

setwd("/IMSB/ra/ibludau/clusterMapper/clusterMapper")
library(devtools)
load_all()
fullData <- copy(data)
fullDataFiltered <- filterData(fullData)
#use_data(fullData,overwrite = TRUE)
use_data(fullDataFiltered,overwrite = TRUE)
use_data(exampleData,overwrite = TRUE)

load_all()
fullDataMapped <- mapClusters(fullDataFiltered)
use_data(fullDataMapped,overwrite = TRUE)


load_all()
data <- exampleData
exampleDataFiltered <- filterData(data,
                                  min_run_coverage = 0.6)
use_data(exampleDataFiltered,overwrite = TRUE)


dataFiltered <- filterData(data,min_run_coverage = 0.6)
complexExample <- subset(dataFiltered,Collapsed_Peptide=="ISSASTPQTSQGRFTAPTSPSTSSPK(3P)")
complexExample <- subset(complexExample,run_id %in% c("0_0","0_1","0_2","0_3","0_4","0_5"))
use_data(complexExample,overwrite = TRUE)


## subset data for code testing

dataMapped <- mapClusters(dataSub, sort_RT=FALSE)
