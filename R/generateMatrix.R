#input <- fread("/IMSB/ra/ibludau/phospho_cluster_mapping/input/combined_clusters_head_all_mscore.txt")
#output <- filterClusters(fullDataMapped)

#' Generate final output matrix
#'
#' @description Merge input data with mapped and filtered results.
#' @param rawInput data.table with original cluster information.
#' @param processedData data.table with mapped and filtered clutsre information.
#' @return data.table Same as input table but including extra columns with
#' new_cluster information.
#' @import data.table
#' @export
#'
#' @examples
#' data <- exampleData
#' filteredData <- filterData(data,
#'                            min_run_coverage = 0.6)
#' mappedData <- mapClusters(filteredData)
#' mappedDataFiltered <- filterClusters(mappedData,
#'                                      min_main_cluster_freq = 0.6)
#' finalData <- generateMatrix(rawInput = data,
#'                             processedData = mappedDataFiltered)
#'
generateMatrix <- function(rawInput, processedData) {
  if (! all(class(rawInput) == c("data.table","data.frame"))) {
    stop("The provided input table \"rawInput\" is not a data.table.")
  }
  if (nrow(rawInput) == 0) {
    stop("The provided input table \"rawInput\" is empty.")
  }
  if (! all(c("transition_group_id","run_id","RT","Collapsed_Peptide","Cluster") %in% names(rawInput))) {
    stop("The provided input table \"rawInput\" does not contain all required columns.")
  }
  if (! all(class(processedData) == c("data.table","data.frame"))) {
    stop("The provided input table \"processedData\" is not a data.table.")
  }
  if (nrow(processedData) == 0) {
    stop("The provided input table \"processedData\" is empty.")
  }
  if (! all(c("transition_group_id","run_id","RT","Collapsed_Peptide","Cluster","new_cluster") %in% names(processedData))) {
    stop("The provided input table \"processedData\" does not contain all required columns.")
  }
  # rename duplicate column names in processed data
  if( length(which(names(processedData) == "Collapsed_Peptide")) != 1) {
    rename_idx <- which(names(processedData) == "Collapsed_Peptide")[2]
    names(processedData)[rename_idx] = "Collapsed_Peptide_2"
  }
  if( length(which(names(processedData) == "run_id")) != 1) {
    rename_idx <- which(names(processedData) == "run_id")[2]
    names(processedData)[rename_idx] = "run_id_2"
  }
  # rename duplicate column names in rawInput
  if( length(which(names(rawInput) == "RT")) != 1) {
    rename_idx <- which(names(rawInput) == "RT")[2]
    names(rawInput)[rename_idx] = "clusterRT"
  }
  # create merged output file
  combi <- merge(rawInput,processedData,all.x=F,all.y=T,by=c("transition_group_id","run_id","RT","Collapsed_Peptide","Cluster"))
  return(combi)
}
