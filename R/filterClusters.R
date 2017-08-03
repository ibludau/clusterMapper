#' Mapped cluster filter
#'
#' @param table data.table with original cluster information
#' @param min_main_cluster_freq Numeric between 0 and 1. Minimum frequency of a
#' transition group id matched to the same cluster. Default is 0.6.
#' @return data.table Same as input table but including extra columns about
#' run coverage information and filtered by \code{min_main_cluster_freq}.
#' @import data.table
#' @export
filterClusters <- function(table, min_main_cluster_freq = 0.6){
  features <- copy(table)
  # number of runs per transition group
  features[, run_count_per_tg := uniqueN(run_id), by = c("transition_group_id")]
  # main cluster per transition group
  getMostFrequentCluster <- function(X){as.numeric(names(sort(table(X),decreasing=TRUE))[1])}
  features[ , `:=` (main_cluster_of_tg = getMostFrequentCluster(new_cluster)) , by = c("Collapsed_Peptide","transition_group_id") ]
  # determine main cluster prevalance
  featuresFiltered <- subset(features, new_cluster == main_cluster_of_tg)
  # determine frequency of unique cluster assignment per transition group
  featuresFiltered[, run_count_per_tg_main_cluster := uniqueN(run_id), by = c("transition_group_id","new_cluster")]
  featuresFiltered[, main_cluster_tg_freq := run_count_per_tg_main_cluster/run_count_per_tg]
  featuresFiltered <- subset(featuresFiltered, main_cluster_tg_freq >= min_main_cluster_freq)
  # remove unassigned clusters
  featuresFiltered <- subset(featuresFiltered, new_cluster != 0)
  return(featuresFiltered)
}
