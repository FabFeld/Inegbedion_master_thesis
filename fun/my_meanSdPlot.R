my_meanSdPlot <- function(normMtx){
  sdp1 <- vsn::meanSdPlot(normMtx$norm, ranks = T, plot = F)
  sdp1 <- sdp1$gg + ggtitle("Mean counts - Ranked")+scale_y_continuous(trans='log10')+ylab("sd (log scale)") +
    theme_bw()
  
  sdp2 <- vsn::meanSdPlot(normMtx$ntd, ranks = T, plot = F)
  sdp2 <- sdp2$gg + ggtitle("log2 (mean counts) - Ranked") +
    theme_bw() #+ ylim(0,3)
  
  sdp3 <- vsn::meanSdPlot(normMtx$nrl.bl, ranks = T, plot = F)
  sdp3 <- sdp3$gg + ggtitle("rlog (mean counts) - Ranked") +
    theme_bw() #+ ylim(0,3)
  
  sdp1b <- vsn::meanSdPlot(normMtx$norm, ranks = F, plot = F) #original scale with ranks=F
  sdp1b <- sdp1b$gg + ggtitle("Mean counts") +
    theme_bw()
  
  sdp2b <- vsn::meanSdPlot(normMtx$ntd, ranks = F, plot = F) #original scale with ranks=F
  sdp2b <- sdp2b$gg + ggtitle("log2 (mean counts)") +
    theme_bw() #+ ylim(0,3)
  
  sdp3b <- vsn::meanSdPlot(normMtx$nrl.bl, ranks = F, plot = F) #original scale with ranks=F
  sdp3b <- sdp3b$gg + ggtitle("rlog (mean counts)") + 
    theme_bw()#+ ylim(0,3)
  
  print(ggarrange(sdp1, sdp2, sdp3, sdp1b, sdp2b, sdp3b,
                  labels = c("A1", "B1", "C1","A2", "B2", "C2"),
                  ncol = 3, nrow =2))
}