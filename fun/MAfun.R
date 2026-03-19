MAfun <- function(res, title="", topN=0, Symbol=F, LFcut= quantile(abs(res$log2FoldChange),Qt)){
  
  x = res
  if(!any(colnames(x)=="baseMean" | colnames(x)=="log2FoldChange")){
    x$mean.ProtAbund <- 2^(x$mean.ProtAbund)
    colnames(x)[grep("mean.ProtAbund", colnames(x))] <- "baseMean"
    colnames(x)[grep("log2FC", colnames(x))] <- "log2FoldChange"
    datType = "prot"
  } else {datType = "rna"}
  
  my_theme <- theme(plot.title = element_text(size = 10, face = 'plain'), #line = element_line(size = .1),
                    axis.text = element_text(size = 10),
                    axis.title = element_text(size = 10),
                    legend.text = element_text(size = 10),
                    text = element_text(size = 7, face = 'plain'))
  
  ma <- ggpubr::ggmaplot(x, fdr = .05, fc = 2^(LFcut), size = 1, top = topN, legend = "bottom",
                         select.top.method = 'fc', # fc or padj
                         palette = c("#B31B21", "#1465AC", "darkgray"),
                         genenames = if(Symbol == F){NULL}else{as.vector(res$SYMBOL)},
                         xlab = if(datType == "prot"){"Mean intensity"}else{bquote(~log[2]~ "(mean expression)")},
                         #font.legend = "bold", font.main = "bold", font.label = c("bold", 11), label.rectangle = F
                         ylab = bquote(~log[2]~ "(fold change)")) +
    ggtitle(paste0(title)) + theme_bw() + theme(legend.position = "bottom")
  
  if(LFcut != 0){
    ma <- ma +
      geom_text(aes(x = max(log2(x$baseMean)*.99), y = LFcut*1.5, label = round( LFcut,2)), size = 3.5) + 
      geom_text(aes(x = max(log2(x$baseMean)*.99), y =-LFcut*1.5, label = round(-LFcut,2)), size = 3.5)
  }
  return(ma + my_theme)
}