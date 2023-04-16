#!/bin/R
####### Author: Jin Woo Oh  ############
###### For CRISPR WG paper  ############
###### Last updated: 04/16/2023#########
########################################

library(ggplot2)
library("cowplot")


###############################################  Fig 3(a) violin plots  ########################################################3
x20=read.table(GATA_HCRFF_x20.log2FC_MN1.bedgraph.peak_overlap', sep="")
x50=read.table(GATA_HCRFF_x50.log2FC_MN1.bedgraph.peak_overlap', sep="")
x100=read.table(GATA_HCRFF_x100.log2FC_MN1.bedgraph.peak_overlap', sep="")
x200=read.table(GATA_HCRFF_x200.log2FC_MN1.bedgraph.peak_overlap', sep="")


x20_T = x20[x20$V5=="TRUE",]$V4
x20_F = x20[x20$V5=="FALSE",]$V4
x50_T = x50[x50$V5=="TRUE",]$V4
x50_F = x50[x50$V5=="FALSE",]$V4
x100_T = x100[x100$V5=="TRUE",]$V4
x100_F = x100[x100$V5=="FALSE",]$V4
x200_T = x200[x200$V5=="TRUE",]$V4
x200_F = x200[x200$V5=="FALSE",]$V4


t20 = t.test(x20_F, x20_T)$statistic
t50 = t.test(x50_F, x50_T)$statistic
t100 = t.test(x100_F, x100_T)$statistic
t200 = t.test(x200_F, x200_T)$statistic


# create a dataset

n1 = c(rep("20x", length(x20_F)), rep("20x", length(x20_T)), rep("50x", length(x50_F)), rep("50x", length(x50_T)), rep("100x", length(x100_F)), rep("100x", length(x100_T)), rep("200x", length(x200_F)), rep("200x", length(x200_T)))
n1 = factor(n1, levels = c("20x", "50x","100x","200x"))

n2 = c(rep("outside GATA1 CRE", length(x20_F)), rep("within GATA1 CRE", length(x20_T)), rep("outside GATA1 CRE", length(x50_F)), rep("within GATA1 CRE", length(x50_T)), rep("outside GATA1 CRE", length(x100_F)), rep("within GATA1 CRE", length(x100_T)), rep("outside GATA1 CRE", length(x200_F)), rep("within GATA1 CRE", length(x200_T)))

cell_cov = c(rep(20, length(x20_F) + length(x20_T)), rep(50, length(x20_F) + length(x20_T)), rep(100, length(x20_F) + length(x20_T)), rep(200, length(x200_F) + length(x200_T)))

within_hit = c(rep(FALSE, length(x20_F)), rep(TRUE, length(x20_T)),rep(FALSE, length(x50_F)),rep(TRUE, length(x50_T)),rep(FALSE, length(x100_F)),rep(TRUE, length(x100_T)),rep(FALSE, length(x200_F)),rep(TRUE, length(x200_T)))


df <- data.frame(
                   value=c(x20_F, x20_T, x50_F, x50_T, x100_F, x100_T, x200_F, x200_T),
                   name1 = n1,
                   name2= n2,
                   within_hit = within_hit,
                   cell_cov = cell_cov
                   )

vplot <- ggplot(data = df, aes(x = name1, y=value, fill=name2)) + theme_classic()+ geom_violin(alpha=0.5,position="identity") + scale_colour_manual(values = c("red","blue")) + xlab("Cell Coverage") +  ylab("guide log2FC effect size") + theme(legend.position = "top", legend.title = element_blank(),text=element_text(size=8), axis.text=element_text(size=8), legend.text=element_text(size=8))

