# Vulcano FUN
VulcFun <- function(res, title="", topN=0, Symbol=F, LFcut= quantile(abs(res$log2FoldChange),Qt)){
  if(LFcut == 0){LFcut = 0.0000001}
  
  if(Symbol == T){
    select <- res[order(res$padj),"SYMBOL"]
  }else{
    select <- rownames(res[order(res$padj),])
  }
  
  my_theme <- theme(plot.title = element_text(size = 10, face = 'plain'), line = element_line(size = .6),
                    axis.text = element_text(size = 10),
                    axis.title = element_text(size = 10),
                    legend.text = element_text(size = 8),
                    text = element_text(size = 8, face = 'plain'))
  
  vu <- EnhancedVolcano::EnhancedVolcano(res, x = "log2FoldChange", y = "padj",
                                         title = paste0(title," [LFcut:",round(LFcut,2),"]"), subtitle = NULL,
                                         lab = if(Symbol == F){rownames(res)}else{res$SYMBOL},
                                         selectLab = if(topN == 0){select[NULL]}else{select[1:topN]}, # select topgenes based on padj value
                                         legendLabels = c('NS','Log2 FC','padj','padj & Log2 FC'),
                                         xlab = bquote(~log[2]~ "(fold change)"), ylab = bquote(~-log[10]~"("~italic(padj)~")"),
                                         FCcutoff = LFcut,   #default 1
                                         pCutoff = .05,      #default 10e-6
                                         labSize = 3.0,
                                         pointSize = 1, #default 0.8
                                         col = c("grey30", "grey30", "royalblue", "red2"),
                                         #shape = c(1, 0, 17, 19),   #default 19, http://sape.inf.usi.ch/quick-reference/ggplot2/shape for details
                                         colAlpha = 0.4, #default 0.5
                                         hline = c(0.01, 0.001), hlineCol = c('grey40','grey55'), hlineType = 'dotted', hlineWidth = 0.6,
                                         gridlines.major = T, gridlines.minor = T, drawConnectors = F,
                                         #widthConnectors = 0.2, colConnectors = 'grey15'
  ) + xlim(-max(abs(res$log2FoldChange))*1.05, max(abs(res$log2FoldChange)*1.05)) +
    ylim(0, max(-log10(res$padj))*1.05) + theme_bw() + theme(legend.position = "bottom")
  
  return(vu + my_theme)
  
}