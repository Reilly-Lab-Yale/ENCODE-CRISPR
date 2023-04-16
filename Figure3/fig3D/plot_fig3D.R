
library(ggplot2)
dat = read.table('/mnt/t/data0/joh27/projects/ENCODE/FCC/crispr_group/final_ENCODE_formatted_data/log2FC/HCRFF_sorting_complexity/seqdepth_cellcov_coanalysis/all_cellcov_seqdepth_bootstrapsample.aurpc', header=TRUE)
dat = dat[dat$seqdepth>1,]
# cellcov seqdepth        sample_id       auprc
#100     1000    10      0.8102106596203971
#100     1000    1       0.8102487739487267
t1 = aggregate(dat$auprc, by = list(dat$cellcov, dat$seqdepth), mean)
t2 = aggregate(dat$auprc, by = list(dat$cellcov, dat$seqdepth), sd)
t3 = aggregate(dat$auprc, by = list(dat$cellcov, dat$seqdepth), min)
t4 = aggregate(dat$auprc, by = list(dat$cellcov, dat$seqdepth), max)

df = data.frame(cellcov = t1$Group.1, seqdepth = t1$Group.2, mean_auprc = t1$x, sd_auprc = t2$x, min_auprc = t3$x, max_auprc = t4$x)
cols = c("20" = rgb(1,0,0, maxColorValue=1), "50" = rgb(0.66,0,0.33, maxColorValue=1),
                                "100" = rgb(0.33,0,0.66, maxColorValue=1), "200" = rgb(0.33,0,1, maxColorValue=1))
linetypes = c("20" = "twodash", "50" = "solid", "100" = "solid", "200" = "solid")
ccsd_plot <-  ggplot(df, aes(x = seqdepth, y = mean_auprc , color = as.factor(cellcov), linetype=as.factor(cellcov))) +
  geom_point() + geom_line()  + geom_errorbar(aes(ymin = mean_auprc - 2.576*sd_auprc/sqrt(10), ymax = mean_auprc + 2.576*sd_auprc/sqrt(10), width=20)) +
  xlab("Simulated sequencing depth") + ylab("AUPRC") + theme_classic() +
  theme(text=element_text(size=8), axis.text=element_text(size=8), legend.position = "none",legend.title = element_blank()) +
  labs(color  = "cellcov") +
  scale_color_manual(values = cols)  + xlim(0,1000) + #scale_x_continuous(breaks = round(seq(0, 1000, by = 100),1) , limits = c(0,1000))
  scale_linetype_manual(values = linetypes)
ccsd_plot

