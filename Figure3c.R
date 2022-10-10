
# load packages; install first if not
#install.packages("tidyverse")
#library(dplyr)
library(tidyverse)

# NDD
mydat <- read.table("data_for_Figure3c.txt", header = T)
plot(mydat$dnMIS ~ mydat$dnLGD, main = "NDD DNMs", xlab = "# dnLGD", ylab = "# dnMIS", pch = 20)
segments(0,0, 100,100, lty = 3, lwd = 2, col = "black")
# change ARID1B dnLGD from 118 to 98 for easy plot

# downstream editing was done in the Adobe Illustrator