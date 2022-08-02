#通路圈图
library(statnet)
library(circlize)
library(dplyr)
library(tidyr)
library(data.table)
library(tidyverse)
library(lubridate)
my.data = read.csv("data2_w_w.csv",row.names=1)
rownames(my.data)
colnames(my.data)
grid.col = NULL
grid.col[c("C5AR1",	"CD4","CSF1R","GPR137B","PTPN2","SLAMF8","SYK","TREM2","VSIG4")] = c("red", "yellow","green", "blue","indianred3","indianred4",
"ivory","ivory1","ivory2")
grid.col[c("P1", "P2", "P3", "P4","P5")] = c("red", "yellow","green", "blue","blue")
circos.par(gap.degree = c(rep(2, nrow(my.data)-1), 10, rep(2, ncol(my.data)-1), 10),
           start.degree = 180)
chordDiagram(my.data,
             directional = TRUE,
             diffHeight = 0.06,
             grid.col = grid.col,link.auto = TRUE,
             transparency = 0.5)
circos.clear()