mapClusters <- function(table, sort_RT=FALSE){
  features <- copy(table)
  ## count number of clusters per run and collapsed peptide
  features[, cluster_count_per_run := uniqueN(Cluster), by = c("run_id","Collapsed_Peptide")]
  ## count total number of clusters across all runs
  features[, cluster_count_total := max(cluster_count_per_run), by = c("Collapsed_Peptide")]
  ## if only one cluster across all runs, new cluster is 1
  features[ , new_cluster := ifelse(cluster_count_total == 1, 1, 0)]
  featuresDone <- subset(features, new_cluster != 0)

  ###############################################################
  ## Continue with data requiring matching ######################
  ###############################################################
  ## subset features to contain only cluster that still need to be matched
  featuresSubset <- subset(features, new_cluster == 0)
  ## tg count per cluster, run and collapsed peptide
  featuresSubset[, tg_count_per_cluster := uniqueN(transition_group_id), by = c("run_id","Collapsed_Peptide","Cluster")]
  ## tg count per run and collapsed peptide
  featuresSubset[, tg_count_per_peptide := uniqueN(transition_group_id), by = c("run_id","Collapsed_Peptide")]

  ###############################################################
  ## Determine reference run for each collapsed peptide #########
  ###############################################################
  ## select one of the remaining runs as reference
  featuresSubsetTop <- copy(featuresSubset)
  featuresSubsetTop <- unique(featuresSubsetTop, by=c("run_id","Collapsed_Peptide"))
  ## select most frequent cluster number
  getMostFrequentCluster <- function(X){as.numeric(names(sort(table(X),decreasing=TRUE))[1])}
  featuresSubsetTop[ , `:=` (freqCluster = getMostFrequentCluster(cluster_count_per_run)) , by = c("Collapsed_Peptide") ]
  ## Select reference run:
  ## most frq cluster count
  featuresSubsetTopSub <- subset(featuresSubsetTop, cluster_count_per_run == freqCluster)
  ## most transition groups
  featuresSubsetTopSub[ , maxTgCount := max(tg_count_per_peptide) , by = .(Collapsed_Peptide) ]
  featuresSubsetTopSub <- subset(featuresSubsetTopSub, tg_count_per_peptide == maxTgCount)
  ## select one of the remaining runs as reference
  featuresSubsetTopSub <- unique(featuresSubsetTopSub, by=c("Collapsed_Peptide"))
  ## This table now contains one reference run per collapsed peptide

  ###############################################################
  ## Generate RT sorted cluster names for each reference run ####
  ###############################################################
  ## Generate a full reference table with a cluster assigned to each transition group base don the RT dimension
  featuresSubset[, unique_key := paste(Collapsed_Peptide,run_id,sep="_")]
  featuresSubset[, meanClusterRT := mean(RT), by=c("unique_key","Cluster")]
  featuresSubset[, clusterRank := frank(meanClusterRT,ties.method="dense"), by=c("unique_key")]

  featuresSubsetTopSub[, unique_key := paste(Collapsed_Peptide,run_id,sep="_")]
  referenceData <- copy(featuresSubset)
  referenceData <- subset(referenceData, unique_key %in% featuresSubsetTopSub$unique_key)
  if (sort_RT) {
    referenceData[, new_cluster := clusterRank]
  } else {
  ## Using clusters sorted by RT might be better, but this wasn't used for original data analysis
  ## Use sort_RT=FALSE to reproduce previous results
    #referenceData[, new_cluster := as.numeric(gsub("cluster","",Cluster))+1]
    referenceData[, new_cluster := clusterRank]
    featuresSubset[, clusterRank := as.numeric(gsub("cluster","",Cluster))+1]
  }

  ###############################################################
  ## Map each test run to the reference #########################
  ###############################################################
  x <- featuresSubset[, mapToReference(.SD,reference = referenceData, sort_RT = FALSE), by=.(Collapsed_Peptide,run_id), .SDcols=c(names(featuresSubset))]
  out <- rbind(featuresDone,x, fill=TRUE)
  return(out)
}

mapToReference <- function(test, reference, sort_RT){
  testSub <- copy(test)
  referenceSub <- copy(reference)
  referenceSub <- subset(referenceSub, Collapsed_Peptide==unique(testSub$Collapsed_Peptide))
  if ((all(referenceSub$cluster_count_per_run == 1)) & (all(testSub$cluster_count_per_run == 1))) {
    testSub[, new_cluster := 1]
    return(testSub)
  } else {
    clusterNum_ref=unique(referenceSub$cluster_count_per_run)
    clusterNum_test=unique(testSub$cluster_count_per_run)
    vote=matrix(data=0,nrow=clusterNum_ref[1],ncol=clusterNum_test[1])
    for (l in 1:clusterNum_ref[1]) { # go through all clusters in reference
      #ref_ids = subset(referenceSub,new_cluster == l)$transition_group_id
      cluster_ref <- unique(referenceSub$new_cluster)
      ref_ids = subset(referenceSub,clusterRank == cluster_ref[l])$transition_group_id
      if(length(ref_ids) < 2) { ################# @TODO think about this!!!!!!!!!!
        next
      }
      #idx_ref_cluster <- which(data$Cluster == ref_cluster[l])
      #cluster_new[Reduce(intersect, list(idx_cp,idx_ref_run,idx_ref_cluster))] = l
      for (m in 1:clusterNum_test[1]) {
        cluster_test <- unique(testSub$clusterRank)
        test_ids = subset(testSub, clusterRank == cluster_test[m])$transition_group_id
        #idx_test_cluster <- which(data$Cluster == test_cluster[m])
        if(length(test_ids) < 2) {
          next
        }
        for (n in 1:length(test_ids)){ # go over all treansition_group_ids in test run
          if (test_ids[n] %in% ref_ids){
            vote[l,m] = vote[l,m]+1
            #print(paste0(vote[l,m], " and ",length(test_ids), " and ",test_ids[n]))
          }
          if (vote[l,m] > length(test_ids)/2) {
            testSub$new_cluster[which(testSub$clusterRank == cluster_test[m])] = cluster_ref[l]
            #if (sort_RT) {
            #  testSub$new_cluster[which(testSub$clusterRank == m)] = l
            #} else {
            #  testSub$new_cluster[which(testSub$clusterRank == m)] = referenceSub$new_cluster[which(referenceSub$new_cluster == l)][1]
            #}
            #print(paste0("break ",l," ",m))
            break
          }
        }
      }
    }
    return(testSub)
  }
}
