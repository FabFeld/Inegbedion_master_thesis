## Log2FC distr.
lfcDistPlot <- function(res, q=Qt, LFcut=quantile(abs(res$log2FoldChange),q), 
                        title="", Ylim = c(0,2000), xlab="Log2(fc)") {
  m = max(abs(res$log2FoldChange))
  if(m > 5){m = 5}
  hist(res$log2FoldChange, breaks= 500, main= paste(substance,title), col= "gray50",
       border= "gray50", xlab=xlab, xlim= c(-m, m), ylim = Ylim)
  abline(v= c(LFcut,-LFcut), col= "dodgerblue", lty= 2, lwd= 1)
  legend(x="topright", text.col = "dodgerblue", bty = "n",
         legend = paste0("LFcut: +/-",round(LFcut, digits = 2),
                         "\n",(1-q)*100,"% Quantile"))
}