myDismtx <- function(mtx, method="euclidean", top=500, title="") {
  if(nrow(mtx) < top){top <- nrow(mtx)}
  message(paste0("Plotting Sample Dist for Var.Top:\t\t",top))
  
  # Subsetting the top genes with max variance
  Xvar <- apply(mtx,1,var)
  X <- t(mtx[names(sort(Xvar, decreasing = T)[1:top]),])
  
  # transpose input, calculate sample distance and create a distance matrix
  sampleDist <- dist(X, method = method) #must be one of "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
  distMtx <- as.matrix(sampleDist)
  #rownames(distMtx) <- paste(coldata$Condition, coldata$Tank, sep="-")
  
  # create annotation object for heatmap
  ann_col <- subset(coldata, select = c("Condition"))#,"Substance"))
  ann_row <- subset(coldata, select = c("Tank"))
  #rownames(ann_row) <- paste(coldata$Mixture, coldata$BioReplicate, sep="-")
  
  colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(15)
  pheatmap::pheatmap(distMtx, angle_col = "90", display_numbers = T,
                     #treeheight_col = 40, #default 50
                     #fontsize = 12, #default 10
                     #fontsize_number = 0.7*12,
                     #cellwidth = 26,
                     #cellheight = 26,
                     drop_levels = T,
                     number_format = if(method == "manhattan"){"%1.0f"}else{"%.1f"},
                     main = paste(title,"- Dist:",method,"- topVar:",top),
                     clustering_distance_rows = sampleDist,
                     clustering_distance_cols = sampleDist,
                     annotation_col = ann_col, # Condition!
                     annotation_row = ann_row, # Tank!
                     annotation_colors = ann_colors,
                     col = colors)
}