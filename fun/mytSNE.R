mytSNE <- function(mtx, top = 500, title = "", IDlabs = F) {
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
    stop("`top` must be at least 2 to compute t-SNE.")
  }
  message(paste0("Plotting t-SNE for ", title, " Var.Top:\t\t\t", top))
  
  # Subsetting the top genes with max variance
  Xvar <- apply(mtx,1,var)
  X <- t(mtx[names(sort(Xvar, decreasing = T)[1:top]),])
  
  # perform Rtsne calc.
  if (nrow(X) < 3) {
    stop("Need at least 3 samples to compute t-SNE.")
  }
  perplex = (nrow(X) - 1) / 3
  set.seed(42) #to make results reproducable! 
  tsne <- Rtsne::Rtsne(X, dims=2, initial_dims=nrow(X), check_duplicates=F, num_threads=2,
                       pca = TRUE, perplexity = perplex, #(should not be bigger than 3 * perplexity < nrow(X) - 1, see details for interpretation)
                       theta = 0.0, #Speed/accuracy trade-off (increase for less accuracy), set to 0.0 for exact TSNE (default: 0.5)
                       partial_pca = F, #(requires the irlba package). This is faster for large input matrices (default: FALSE)
                       max_iter = 10000, #number of iterations
                       is_distance = FALSE, pca_center = TRUE, pca_scale = FALSE, verbose = F,
                       normalize = F) #Default True; Set to F as DESeq's RLE norm was performed prior!
  row.names(tsne$Y) <- row.names(X)
  colnames(tsne$Y) <- c("tSNE_1","tSNE_2")
  
  # ggplot
  ggtsne <- merge(coldata, tsne$Y, by=0)
  shape_var <- "Tank"
  if (!"Tank" %in% names(ggtsne) || all(is.na(ggtsne$Tank))) {
    shape_var <- "Condition"
  }
  if (!exists("ann_colors", inherits = TRUE) || is.null(ann_colors$Condition)) {
    stop("`ann_colors$Condition` not found. Provide colors for `Condition`.")
  }
  g = ggplot(ggtsne) +
    geom_point(aes(x = tSNE_1, y = tSNE_2, color = Condition, shape = .data[[shape_var]]), size = 3, alpha = .7) +
    ggtitle(paste(title, "; Top:", top)) + xlab("t-SNE 1") + ylab("t-SNE 2") +
    scale_colour_manual(values = ann_colors$Condition) + theme_bw() + theme(aspect.ratio = 1)
  if(IDlabs == T){
    return(g + ggrepel::geom_text_repel(aes(x=tSNE_1, y=tSNE_2, label=Row.names), nudge_x = .5, nudge_y = .5))
  } else {
    return(g)
  }
}