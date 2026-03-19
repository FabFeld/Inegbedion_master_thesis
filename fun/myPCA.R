myPCA  <- function(mtx, pcaM = "svd", top = 500, title = "", IDlabs = F) {
  if (!is.matrix(mtx) && !is.data.frame(mtx)) {
    stop("`mtx` must be a matrix or data.frame.")
  }
  if (!exists("coldata", inherits = TRUE)) {
    stop("`coldata` not found. Provide a `coldata` data.frame in scope.")
  }
  if (!is.data.frame(coldata)) {
    stop("`coldata` must be a data.frame.")
  }
  if (!"Condition" %in% names(coldata)) {
    stop("`coldata` must contain a `Condition` column.")
  }
  if (nrow(mtx) < top) {
    top <- nrow(mtx)
  }
  if (top < 2) {
    stop("`top` must be at least 2 to compute PCA.")
  }
  message(paste0("Plotting PCA for ", title, " Var.Top:\t\t\t", top))
  
  # Subsetting the top genes with max variance
  Xvar <- apply(mtx,1,var)
  X <- t(mtx[names(sort(Xvar, decreasing = T)[1:top]),])
  pc <- pcaMethods::pca(X, method = pcaM, center = T, nPcs=2)
  pcaDf <- merge(coldata, pcaMethods::scores(pc), by=0)
  
  shape_var <- "Tank"
  if (!"Tank" %in% names(pcaDf) || all(is.na(pcaDf$Tank))) {
    shape_var <- "Condition"
  }
  if (!exists("ann_colors", inherits = TRUE) || is.null(ann_colors$Condition)) {
    stop("`ann_colors$Condition` not found. Provide colors for `Condition`.")
  }
  
  g = ggplot(pcaDf, aes(PC1, PC2, colour = Condition, shape = .data[[shape_var]])) +
    geom_point(size = 3, alpha = .65) + scale_colour_manual(values = ann_colors$Condition) +
    ggtitle(paste0(title," - exVar: ",round(pc@R2cum[2]*100,1),"%; ",pcaM,"; Top:",pc@nVar)) +
    xlab(paste0("PC1: ",round((pc@R2)[1]*100,1),"% variance")) +
    ylab(paste0("PC2: ",round((pc@R2)[2]*100,1),"% variance")) + #stat_ellipse() +
    theme_bw() + theme(aspect.ratio = 1, legend.position = "none") +
    theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12))
  if(IDlabs == T){
    return(g + ggrepel::geom_text_repel(aes(x=PC1, y=PC2, label=Row.names), nudge_x = .5, nudge_y = .5))
  } else {
    return(g)
  }
}