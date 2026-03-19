DE.Cor.Venn = function(df.x, df.y, pcut = .05, typ.x = "DEG in HE", typ.y = "DEG in LE", title = "",
                       LFcut.x = quantile(abs(df.x$log2FC),Qt),
                       LFcut.y = quantile(abs(df.y$log2FC),Qt),
                       COL = c("red3","deepskyblue1","blue4"),
                       lfsSelect = F, lfsSet.x = NULL, lfsSet.y = NULL){
  # subset df.x/.y for plotting based on pcut and LFcut settings
  x = subset(df.x[df.x$padj <= pcut, ], abs(log2FC) >= LFcut.x)[,"Gene"] #%>% droplevels(.)
  y = subset(df.y[df.y$padj <= pcut, ], abs(log2FC) >= LFcut.y)[,"Gene"] #%>% droplevels(.)
  
  # Provide custom gene set selection based on lfs cut
  if(lfsSelect == T){
    message("LF shrunk DEG selection provided")
    x = lfsSet.x
    y = lfsSet.y
  }
  
  ## PLOT ONLY IF x & y contains data!!!
  if( length(x) + length(y) > 0){
    
    # get common DEGs
    ls = list(DEG.x = unique(x), DEG.y = unique(y))
    names(ls) <- c(typ.x, typ.y)
    int = intersect(ls[[1]], ls[[2]])
    uni = union(ls[[1]], ls[[2]]) # union set of x & y
    
    # write df for gg corrplot
    X = df.x[df.x$Gene %in% uni,c("Gene","log2FC")] #%>% droplevels(.)
    Y = df.y[df.y$Gene %in% uni,c("Gene","log2FC")] #%>% droplevels(.)
    stopifnot(nrow(X) == nrow(Y))
    df = merge(X,Y, by = "Gene")[,-4]
    
    # append Type info to df
    df[df$Gene %in% ls[[1]], "Type"] <- typ.x
    df[df$Gene %in% ls[[2]], "Type"] <- typ.y
    df[df$Gene %in% int, "Type"] <- "DEG Overlap"
    df$Type <- factor(df$Type, levels = c("DEG Overlap", typ.y, typ.x))
    
    col = setNames(COL,c("DEG Overlap", typ.y, typ.x))
    gg = list()
    
    ## Venn plot ##
    venn = eulerr::euler(ls[c(2,1)], shape = "ellipse")
    s = round(venn$stress,3)
    e = round(venn$diagError,3)
    tmp = sort(unlist(lapply(ls, length)))
    pct = round(length(int)/tmp[1]*100,1)
    gg[["VE"]] = plot(venn,
                      fills = list(fill = col[c(2,3,1)], alpha = .6), #labels = list(col = "black", font = 4),
                      legend = list(col = "black", font = 4),
                      main = paste0("Stress: ",s,"; Diag.Er: ",e,"\nShared ",pct,"%"), #main = paste0(fn,": sign. terms [padj < ",pcut,"]"),
                      quantities = TRUE, shape = "ellipse", lty = 0)
    
    ## Corr plot ##
    # cor test
    cdf <- df[which(df$Type %in% 'DEG Overlap'),]
    if(nrow(cdf) > 2){
      res <- cor.test(cdf$log2FC.x, cdf$log2FC.y, method = "pearson", alternative = "greater")
      p <- signif(res$p.value, digits = 2)
      t <- round(res$statistic,1)
      d <- round(res$parameter,1)
      c <- round(res$estimate,2)
      if(p <= .001) {i = "***"
      }else if(p <= .01) {i = "**"
      }else if(p <= .05) {i = "*"
      }else {i = "NS"}
    } else {
      warning("Not enough common observations between HE & LE to run cor.test!\n")
      c <- round(cor(cdf$log2FC.x, cdf$log2FC.y, method = "pearson"),2)
      p <- t <- d <- i <- "NA"
    }
    cA <- round(cor(df$log2FC.x,df$log2FC.y),2)
    
    my_theme <- theme(axis.text = element_text(size = 12),
                      axis.title = element_text(size = 13))
    # cor plot
    gg[["COR"]] = ggplot(df, aes(x=log2FC.x, y=log2FC.y, color = Type)) +
      geom_point(size = 3, alpha = .4, shape = 16) +
      geom_point(data = df[which(df$Type %in% 'DEG Overlap'),], #2nd layer
                 aes(x=log2FC.x, y=log2FC.y), shape = 16, size = 3, alpha = .5) +
      scale_color_manual(values = col) +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = .4) +
      geom_vline(xintercept = 0, linetype = "dashed", alpha = .4) +
      theme_bw() + theme(legend.position = 'bottom') + geom_rug(alpha = .5) +
      labs(title = paste(substance,"- DEG cor:",title), 
           subtitle = paste0(typ.x," vs ", typ.y,": (t=",t," df=",d," p=",p,")",i), 
           x = paste0("Log2FC (",typ.x,")"), 
           y = paste0("Log2FC (",typ.y,")")) +
      annotate("text", label = paste0("Pearson Cor = ",c,"\nR2 = ",round(c^2,2)),
               color = COL[1],
               x = (min(df$log2FC.x)), 
               y = (max(df$log2FC.y)),
               hjust=0, vjust=.52) + 
      annotate("text", label = paste0("Pearson Cor = ",cA,"\nR2 = ",round(cA^2,2)),
               x = (max(df$log2FC.x)), 
               y = (min(df$log2FC.y)),
               hjust=.9, vjust=0) + my_theme
    
    return(gg)
  } else {
    return(NULL)
  }
}