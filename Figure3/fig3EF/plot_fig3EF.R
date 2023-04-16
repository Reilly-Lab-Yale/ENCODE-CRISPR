####### Author: Jin Woo Oh  ############
###### For CRISPR WG paper  ############
###### Last updated: 04/16/2023#########
########################################

################################################### Fig 3(e) #####################################################################
library(ggplot2); library(reshape2)
# GATA1 Bassik
gata1_i_c = read.table('I_Bassik_GATA1.samp_corrs_stat',col.names=c('depth', 'mean', 'std', 'min', 'max'))
gata1_k_c = read.table('K_Bassik_GATA1.samp_corrs_stat', col.names=c('depth', 'mean', 'std', 'min', 'max'))

# Doughty data
gata1_eg_c = read.table('Engreitz_GATA1_growth.samp_corrs_stat', sep=""),col.names=c('depth', 'mean', 'std', 'min', 'max'))
gata1_ef_c = read.table('FF_GATA1.samp_corrs_stat' ,col.names=c('depth', 'mean', 'std', 'min', 'max'))
myc_eg_c = read.table('Engreitz_MYC_growth.samp_corrs_stat', ,col.names=c('depth', 'mean', 'std', 'min', 'max'))

# Reily data
gata1_r_c = read.table(paste(dir,'HCRFF_GATA1.samp_corrs_stat', sep=""), ,col.names=c('depth', 'mean', 'std', 'min', 'max'))
myc_r_c = read.table(paste(dir,'HCRFF_MYC.samp_corrs_stat', sep=""),,col.names=c('depth', 'mean', 'std', 'min', 'max'))

table_vector_c = vector(mode = "list", length = 7)
table_vector_c[[1]] = gata1_i_c
table_vector_c[[2]] = gata1_k_c
table_vector_c[[3]] = gata1_eg_c
table_vector_c[[4]] = gata1_ef_c
table_vector_c[[5]] = myc_eg_c
table_vector_c[[6]] = gata1_r_c
table_vector_c[[7]] = myc_r_c
# normalize by maximum seq depth (x5000)
for(i in 1:7){
        dat = table_vector_c[[i]]
        max_corr = dat[dim(dat)[1],]$max
        table_vector_c[[i]]$min = table_vector_c[[i]]$min / max_corr
        table_vector_c[[i]]$max = table_vector_c[[i]]$max / max_corr
        table_vector_c[[i]]$mean = table_vector_c[[i]]$mean / max_corr
        table_vector_c[[i]]$std = table_vector_c[[i]]$std / max_corr
}

ds = 1
es = .1
cplot = ggplot(data=table_vector_c[[1]], aes(x = depth,y=mean))+
        geom_point(data=table_vector_c[[1]], aes(y=mean, color = method[1], shape = gene[1])) +

        geom_point(data=table_vector_c[[3]], aes(y=mean,  color = method[3], shape = gene[3])) +

        geom_point(data=table_vector_c[[4]], aes(y=mean,  color = method[4], shape = gene[4])) +

        geom_point(data=table_vector_c[[5]], aes(y=mean,  color = method[5], shape = gene[5])) +

        geom_point(data=table_vector_c[[6]], aes(y=mean,  color = method[6], shape = gene[6])) +

        geom_point(data=table_vector_c[[7]], aes(y=mean,  color = method[7], shape = gene[7]))

cplot <- cplot+ theme_classic() + scale_color_manual(values = colors) + scale_shape_manual(values = shapes) + theme(legend.position = "none",legend.title = element_blank(),text=element_text(size=8), axis.text=element_text(size=8)) + ylab("Replicate log2FC corr (relative to x5000)") + xlab("Simulated sequencing depth") + xlim(0,1000) + ylim(0.5,1)









##############################################  Fig 3(f)  #############################################################

# GATA1 Bassik
gata1_i_d = read.table('I_Bassik_GATA1.samp_dropout_stat', col.names=c('depth', 'mean', 'std', 'min', 'max'))
gata1_k_d = read.table('K_Bassik_GATA1.samp_dropout_stat',  col.names=c('depth', 'mean', 'std', 'min', 'max'))
# Doughty data
gata1_eg_d = read.table('Engreitz_GATA1_growth.samp_dropout_stat',col.names=c('depth', 'mean', 'std', 'min', 'max'))
gata1_ef_d = read.table('FF_GATA1.samp_dropout_stat',col.names=c('depth', 'mean', 'std', 'min', 'max'))
myc_eg_d = read.table('Engreitz_MYC_growth.samp_dropout_stat',col.names=c('depth', 'mean', 'std', 'min', 'max'))
# Reily data
gata1_r_d = read.table('HCRFF_GATA1.samp_dropout_stat',col.names=c('depth', 'mean', 'std', 'min', 'max'))
myc_r_d = read.table('HCRFF_MYC.samp_dropout_stat',col.names=c('depth', 'mean', 'std', 'min', 'max'))

table_vector_d = vector(mode = "list", length = 15)
table_vector_d[[1]] = gata1_i_d
table_vector_d[[2]] = gata1_k_d
table_vector_d[[3]] = gata1_eg_d
table_vector_d[[4]] = gata1_ef_d
table_vector_d[[5]] = myc_eg_d
table_vector_d[[6]] = gata1_r_d
table_vector_d[[7]] = myc_r_d

for(i in 1:7){
        table_vector_d[[i]]$mean = table_vector_d[[i]]$mean * 100
        table_vector_d[[i]]$std = table_vector_d[[i]]$std * 100
        table_vector_d[[i]]$min = table_vector_d[[i]]$min * 100
        table_vector_d[[i]]$max = table_vector_d[[i]]$max * 100

}


ds = 1
es = .1
dplot = ggplot(data=table_vector_d[[1]], aes(x = depth,y=mean))+
        geom_point(data=table_vector_d[[1]], aes(y=mean, color = method[1], shape = gene[1])) +

        geom_point(data=table_vector_d[[3]], aes(y=mean,  color = method[3], shape = gene[3])) +

        geom_point(data=table_vector_d[[4]], aes(y=mean,  color = method[4], shape = gene[4])) +

        geom_point(data=table_vector_d[[5]], aes(y=mean,  color = method[5], shape = gene[5])) +

        geom_point(data=table_vector_d[[6]], aes(y=mean,  color = method[6], shape = gene[6])) +

        geom_point(data=table_vector_d[[7]], aes(y=mean,  color = method[7], shape = gene[7]))

dplot <- dplot + theme_classic()+ scale_color_manual(values = colors) + scale_shape_manual(values = shapes) + theme(legend.position = "none",legend.title = element_blank() ,text=element_text(size=8), axis.text=element_text(size=8)) + ylim(0,40) + ylab("Dropout rate (%)") + xlab("Simulated sequencing depth") + xlim(0,1000)



