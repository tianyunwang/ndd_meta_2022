# Code to plot Figure2b
# install.packages("UpSetR")
library(UpSetR)

# Dataset
FDR <- c(
  CH=91,
  DNW=142,
  DR=33,
  "CH&DNW"=59,
  "CH&DR"=12,
  "DR&DNW"=41,
  "CH&DR&DNW"=237
)

FWER <- c(
  CH=51,
  DNW=32,
  DR=14,
  "CH&DNW"=10,
  "CH&DR"=9,
  "DR&DNW"=10,
  "CH&DR&DNW"=138
)

upset(fromExpression(FDR), 
      sets.bar.color = "gray",
      matrix.color = "black", 
      main.bar.color = "black",
      nintersects = 40, 
      nsets = 6, 
      order.by = "freq", 
      decreasing = T, 
      mb.ratio = c(0.7, 0.3),
      number.angles = 0, 
      text.scale = 2, 
      point.size = 4, 
      line.size = 1,
      sets.x.label = "# genes",
      mainbar.y.label = "# genes"
)

upset(fromExpression(FWER), 
      sets.bar.color = "gray",
      matrix.color = "black", 
      main.bar.color = "red",
      nintersects = 40, 
      nsets = 6, 
      order.by = "freq", 
      decreasing = T, 
      mb.ratio = c(0.7, 0.3),
      number.angles = 0, 
      text.scale = 2, 
      point.size = 4, 
      line.size = 1,
      sets.x.label = "# genes",
      mainbar.y.label = "# genes"
)

# downstream editing was done in the Adobe Illustrator

