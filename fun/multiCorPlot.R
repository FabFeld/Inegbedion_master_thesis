multiCorPlot = function(mtx, label="") {
  gg.ls = list()
  for(k in names(comb)){
    x = comb[[k]]
    for(i in 1:nrow(x)){
      tmp = x[i,] #IDs to plot with
      df = as.data.frame(mtx[ ,tmp])
      gg.ls[[paste0(k,i)]] <- corplot(df, tmp, Lab = label,
                                        title = gsub("[Ee]xposure","",k))
    }
  }
  gg.ls
}