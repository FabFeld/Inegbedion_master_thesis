# this function takes the raw countmatrix and filtered count matrix as input and
# plots histograms for individual counts and row mean counts


geneCountDistribution <- function(rawCounts, filteredCounts, breaks = 100){
  tmp  = apply(count.matrix,1,mean)
  tmp1 = apply(count.matrix.cl, 1, mean)
  
  h <- hist(as.matrix(log2(tmp+1)), breaks = breaks, plot = F)
  h2 <- hist(as.matrix(log2(count.matrix+1)), breaks = breaks, plot = F)
  
  par(mfrow = c(2, 2))
  
  hist(log2(tmp+1), breaks = breaks, main = "countMtx", 
       ylim = c(0, max(h$counts)),
       xlim = c(0,20),
       xlab = "log2(Row mean count +1)")
  hist(log2(as.matrix(count.matrix+1)), breaks = breaks, main = "countMtx", 
       ylim = c(0, max(h2$counts)),
       xlim = c(0,20),
       xlab = "log2(counts +1)")
  hist(log2(tmp1+1), breaks = breaks, main = "filtered countMtx", 
       ylim = c(0, max(h$counts)),
       xlim = c(0,20),
       xlab = "log2(Row mean count +1)")
  hist(log2(as.matrix(count.matrix.cl+1)), breaks = breaks, main = "filtered countMtx", 
       ylim = c(0, max(h2$counts)),
       xlim = c(0,20),
       xlab = "log2(counts +1)")
  
}


