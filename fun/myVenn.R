myVenn = function(deg, title = "", shape = "ellipse", ...) {
  # Skip plotting when no DEGs in the input
  if(all(sapply(lapply(deg,length),sum) == 0)){
    warning("\nNo DEGs found in any treatment condition. Skipping multi-Venn plotting!")
  } else {
    # Venn fun
    set.seed(42)
    venn <- eulerr::euler(deg, shape = shape, ...)
    s <- round(venn$stress,3)
    e <- round(venn$diagError,3)
    
    # plot
    return(
      plot(venn,
           fills = list(fill = ann_colors$Condition[-1], alpha = .6),
           legend = list(col = "black", font = 4),
           main = paste0(title,"\nStress: ",s," - Diag.Er: ",e),
           quantities = TRUE, shape = shape, lty = 0)
    )
  }
}