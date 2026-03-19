## pval vs padj & qval
pVSpadjPlot <- function(res, pval = 0.05, q = 0.05, p = 0.05, title = "") {
  plot(res$pvalue, res$padj, xlab = "pvalue", ylab = "conversion of pvalue",
       xlim = c(0, 0.25), ylim = c(0, 1), main = paste(substance, title),
       col = "dodgerblue1", pch = 20)
  points(res$pvalue, res$qval, col = "springgreen3", pch = 20)
  abline(h = p, lwd = 1, lty = 2)
  abline(v = pval, lwd = 1, lty = 2)
  xpos <- 0.1
  text(x = 0.25, y = 0.2, pos = 2, offset = 0, col = "black",
       labels = paste("Genes with pvalue <", pval, ":", sum(res$pvalue <= pval, na.rm = TRUE)))
  text(x = 0.25, y = 0.15, pos = 2, offset = 0, col = "springgreen3",
       labels = paste("Genes with Qvalue <", q, ":", sum(res$qval <= q, na.rm = TRUE)))
  text(x = 0.25, y = 0.1, pos = 2, offset = 0, col = "dodgerblue",
       labels = paste("Genes with padj (BH) <", p, ":", sum(res$padj <= p, na.rm = TRUE)))
}
