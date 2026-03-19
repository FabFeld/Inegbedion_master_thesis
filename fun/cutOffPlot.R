# rowsum cutoff vs remaining number of genes in count matrix
cutOffPlot <- function(countMtx, cut) {
  n = ncol(countMtx)
  if(missing(cut)){cut = n}
  if(n < 6) {stop("Number of count matrix columns < 6")}
  X = seq(n-6,n+50,1) 
  X[X == 0] <- 1
  X = sort(union(X, c(1:5)))
  
  rNbr = c() # empty vector to store Nbr of genes per cutoff in
  for(x in X) {
    tmp  = nrow(countMtx[which(rowSums(countMtx) >= x), ])
    rNbr = append(rNbr, tmp)
  }
  df <- data.frame(x = X, y = rNbr)
  plot(df, xlab = "Row sum cut off", ylab = "Nbr of genes",
       main = "Genes in count matrix - Deprecated, do not use")
  abline(v=cut, lwd=1.5, lty=2, col = "firebrick")
  #abline(h=df[cut,2], lwd=1.5, lty=2, col = "blue")
  text(x = cut, y = df[cut,2], pos = 4, offset = 1.5,
       labels = paste("Cutoff:",cut,"=",df[cut,2],"genes"))
}