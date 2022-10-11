###########################
# plot for Figure 3a and 3b
###########################
DNM <- read.delim("data_for_Figure3a_and_3b.txt", header = T)

# plot for Figure 3a - dnLGD frequency
plot(DNM$dnLGD.ASD.pct ~ DNM$dnLGD.DD.pct, main = "dnLGD frequency", 
     xlab = "DD (n=31,052)", ylab = "ASD (n=15,560)", pch = 20, col = as.numeric(DNM$Sig.col))
# change to plot ARID1B dnLGD.DD.pct from 0.34 to 0.29 for easy plot
segments(-0.01,-0.01, 0.14,0.14, lty = 3, lwd = 2, col = "black")

# plot for Figure 3a - dnMIS frequency
plot(DNM$dnMIS.ASD.pct ~ DNM$dnMIS.DD.pct, main = "dnMIS frequency", 
     xlab = "DD (n=31,052)", ylab = "ASD (n=15,560)", pch = 16, col = as.numeric(DNM$Sig.col))
segments(-0.01,-0.01, 0.14,0.14,  lty = 3, lwd = 2, col = "black")

# downstream editing was done in the Adobe Illustrator