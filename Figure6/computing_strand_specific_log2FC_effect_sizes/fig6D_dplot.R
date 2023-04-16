library(ggplot2)
library("cowplot")

colors = c('growth' = 'purple' , 'HCRFF' = 'black')
shapes = c('CRISPRi' = 19, 'CRISPRa' = 17, 'CRISPRd' = 1 , 'CRISPRk' = 15)


dat = read.table('samp_to_mean_log2fc_z_by_type_ravg.txt', header=TRUE, sep="\t")

dat = dat[dat$N_cgb > 10 & dat$N_cp > 10 & dat$N_tgb > 10 & dat$N_tp > 10, ]
#samp    perturbation_type       readout coding_gene_body        N_cgb   template_gene_body      N_tgb   coding_promoter N_cp    template_promoter       N_tp
#A_Bassik        CRISPRa growth  0.2750521583901774      242     0.04626102865049488     242     -0.33874846888534105    78      -0.28666676265103314    57

df <- data.frame(val_coding = dat$coding_gene_body, val_template = dat$template_gene_body, perturbation_type = dat$perturbation_type, readout = dat$readout)

dplot1 <- ggplot(data = df, aes(x = val_template, y = val_coding,  color = readout, shape = perturbation_type)) + theme_classic()+ geom_point() +  xlab("avg. template guide log2FC (Z)") +  ylab("avg. coding log2FC (Z)") + scale_color_manual(values = colors) + scale_shape_manual(values = shapes) + theme(legend.position = "none", text=element_text(size=8), axis.text=element_text(size=8)) + xlim(min(c(dat$template_gene_body, dat$coding_gene_body)), max(c(dat$template_gene_body, dat$coding_gene_body)))  + geom_abline(slope=1, intercept=0)


dplot1
