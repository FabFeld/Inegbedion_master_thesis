# plotting function
corplot = function(df, tmp, Lab="", title="") {
  d  = densCols(df[,tmp[1]], df[,tmp[2]], colramp = colorRampPalette(rev(rainbow(10, end = 4/6))))
  ggplot(df) + geom_point(aes(df[,tmp[1]], df[,tmp[2]], col = d), size = .8, alpha = .4) +
    labs(x=paste(Lab,tmp[1]),y=paste(Lab,tmp[2]),
         subtitle = (paste0(title," [",coldata[tmp[1],"Tank"]," vs ",coldata[tmp[2],"Tank"],"]"))) +
    annotate("text", label = paste0("Pearson: ",round(cor(df[,tmp[1]], df[,tmp[2]]),2)),
             x = (max(df[,tmp[1]])), y = (min(df[,tmp[2]])), hjust=1, vjust=0) + 
    annotate("text", label = paste0("R2: ",round((cor(df[,tmp[1]], df[,tmp[2]]))^2,2)),
             x = (min(df[,tmp[1]])), y = (max(df[,tmp[2]])), hjust=0, vjust=1) +
    scale_color_identity() + coord_equal(ratio=1) + theme_bw()
}