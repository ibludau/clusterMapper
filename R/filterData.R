#' Coverage filter
#'
#' @param table data.table with original cluster information
#' @param min_run_coverage Numeric between 0 and 1. Minimum coverage of a
#' transition group id across all runs. Default is 0.6.
#' @return data.table Same as input table but including extra columns about
#' run coverage information and filtered by \code{min_run_coverage}.
#' @import data.table
#' @export
#'
#' @examples
#' data <- exampleData
#' filteredData <- filterData(data,
#'                            min_run_coverage = 0.6)
filterData <- function(table,
                       min_run_coverage = 0.6){
  if (! all(class(table) == c("data.table","data.frame"))) {
    stop("The provided input table is not a data.table.")
  }
  if (nrow(table) == 0) {
    stop("The provided input table is empty.")
  }
  if (! all(c("transition_group_id","run_id","RT","Collapsed_Peptide","Cluster") %in% names(table))) {
    stop("The provided input table does not contain all required columns.")
  }
  if ((min_run_coverage > 1) | (min_run_coverage < 0)) {
    stop("The provided min_run_coverage value is not valid. Please choose a value between 0 and 1.")
  }
  features <- copy(table)
  features <- subset(features, select = c("transition_group_id","run_id","RT","Collapsed_Peptide","Cluster"))
  features[, run_freq := .N , by = transition_group_id ]
  runCount <- length(unique(features$run_id))
  minRunCount <- ceiling(runCount*min_run_coverage)
  featuresFiltered <- subset(features, run_freq >= minRunCount)
  featuresRemoved <- subset(features, run_freq < minRunCount)
  removedIds <- length(unique(featuresRemoved$transition_group_id))
  totalIds <- length(unique(features$transition_group_id))
  percentageRemoved <- round(100/totalIds*removedIds,digits = 2)
  message(paste0("Removed ",removedIds," out of ",totalIds,
                 " transition_group_ids (",percentageRemoved,
                 "%) based on a min_coverage_filter of ",
                 min_run_coverage,"."))
  return(featuresFiltered)
}

