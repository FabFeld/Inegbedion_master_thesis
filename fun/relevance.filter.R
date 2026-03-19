# relevance filter function for raw count matrix table.
# performs CPM transformation and subsequently filtering and returns 
# non-transformed but filtered count mtx)


relevance.filter = function(mtx, id.ls, minCPM=1){
  
  # Check the number of levels in coldata$Tank (= replicates)
  #num_levels <- length(levels(coldata$Tank))
  num_levels <- length(coldata$Tank)/length(levels(coldata$condition))
  
  # Set the value of perc based on the number replicates
  if (num_levels == 3) {
    perc <- 0.65
  } else if (num_levels >= 4) {
    perc <- 0.75
  }
  mtx_CPM <- my_CPM(mtx) # transform matrix in counts per million
  
  # subset mtx for each exp. condition 
  cond = lapply(id.ls, function(x){ mtx_CPM[,x] })
  
  # check if a certain percentage (perc.) of CPMs is > minCPM
  cond = lapply(cond, function(x){
    X = apply(x, 1, function(i){
      # Count the values > minCPM in i and divide by the length of i
      Z = length(which((i >= minCPM) == T)) / length(i)
      # If Z >= the min percentage criteria return T else F
      if(Z >= perc){ TRUE } else { FALSE }
    })
    names(X) <- row.names(x)
    return(X)
  })
  
  # Combine all in a df
  cond = as.data.frame(rlist::list.cbind(cond))
  
  # retrieve only GeneIDs for which atleast one condition with all (75%) replicates >1 CPM
  keep = names(which(apply(cond,1, function(x){ any(x == T) }) == T))
  
  # Calculate the percentage of kept gene IDs
  n = round((length(keep) / nrow(mtx_CPM)) * 100, 2)
  
  # Print the number of gene IDs that met the relevance criteria
  # and the percentage of initial gene IDs that were kept for downstream processing
  message(length(keep), " out of ", nrow(mtx_CPM), " gene IDs met the relevance criteria.\n",
          n, "% of initial gene IDs were kept for downstream processing.")
  
  # Return the filtered matrix based on the kept gene IDs
  return(mtx[keep,])
}
