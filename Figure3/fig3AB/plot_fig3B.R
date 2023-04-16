#!/bin/R
library(ggplot2)
library("cowplot")

#PRC
prc20 = read.table('precision_recall_x20.out'); colnames(prc20) = c("P", "R")
prc50 = read.table('precision_recall_x50.out');  colnames(prc50) = c("P", "R")
prc100 = read.table('precision_recall_x100.out');  colnames(prc100) = c("P", "R")
prc200 = read.table('precision_recall_x200.out');  colnames(prc200) = c("P", "R")

library(reshape2); library(ggplot2)
prc20.melt<-melt(prc20[,c("P", "R")], id="P")
prc50.melt<-melt(prc50[,c("P", "R")], id="P")
prc100.melt<-melt(prc100[,c("P", "R")], id="P")
prc200.melt<-melt(prc200[,c("P", "R")], id="P")

prc20.melt[, "coverage"] <-"x20"
prc50.melt[, "coverage"] <-"x50"
prc100.melt[, "coverage"] <-"x100"
prc200.melt[, "coverage"] <-"x200"

df <- rbind(prc20.melt, prc50.melt, prc100.melt, prc200.melt)

cols = c("x20" = rgb(1,0,0, maxColorValue=1), "x50" = rgb(0.66,0,0.33, maxColorValue=1),
                                "x100" = rgb(0.33,0,0.66, maxColorValue=1), "x200" = rgb(0.33,0,1, maxColorValue=1))


prc_plot <- ggplot() +
  geom_line(data = df, aes(x = P, y = value, color = coverage, linetype = coverage)) +
  xlab("precision") +  ylab("recall") + theme_classic() +
  theme(text=element_text(size=8), axis.text=element_text(size=8), legend.text=element_text(size=8)) +
  labs(color  = "coverage", linetype = "coverage") +
  scale_color_manual(values = cols) +
  scale_linetype_manual(values = c("x20" = "twodash", "x50" = "solid", "x100" = "solid","x200" = "solid" ))

