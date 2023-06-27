
### Load packages

suppressPackageStartupMessages(library(tidyverse))

`%notin%` <- Negate(`%in%`)

### Load and process perturbation data

#### Bin coordinates with perturbations

allregions <- read.delim("/data/gersbachlab/lrb53/encodeCrisprWG/finalversion/figure1/k562meta/elementReference/combined/binnedgenome.int.k562.elementReference.20230415.merge.bed",
                         header=FALSE, col.names=c("bin.chr","bin.start","bin.end")) %>%
mutate(bin = paste0(bin.chr,":",bin.start,"-",bin.end)) %>% select(bin) %>% distinct()

head(allregions,n=1)
dim(allregions)

#### Bin coordinates with CREs

sigregions <- read.delim("/data/gersbachlab/lrb53/encodeCrisprWG/finalversion/figure1/k562meta/elementQuant/combined/binnedgenome.int.k562.elementQuant.20230415.merge.bed",
                         header=FALSE, col.names=c("bin.chr","bin.start","bin.end")) %>%
mutate(bin = paste0(bin.chr,":",bin.start,"-",bin.end)) %>% select(bin) %>% distinct()

head(sigregions,n=1)
dim(sigregions)

### Load intersection with `bedgraph` files for features

##### All regions

allregion_h3k27ac <- read.delim("/data/gersbachlab/lrb53/encodeCrisprWG/finalversion/figure1/k562meta/outs/comparesignal/binnedgenome.int.k562.elementReference.int.k562.h3k27ac.foldchange.txt",
                             header = FALSE, col.names = c("feature.chr","feature.start","feature.end","signal","bin.chr","bin.start","bin.end")) %>%
mutate(bin = paste0(bin.chr,":",bin.start,"-",bin.end),
       feature_coords = paste0(feature.chr,":",feature.start,"-",feature.end),
       feature = "H3K27ac") %>% 
select(feature, bin, feature_coords, signal) %>% distinct()

head(allregion_h3k27ac,n=1)
dim(allregion_h3k27ac)

allregion_dnase <- read.delim("/data/gersbachlab/lrb53/encodeCrisprWG/finalversion/figure1/k562meta/outs/comparesignal/binnedgenome.int.k562.elementReference.int.k562.dnase.foldchange.txt",
                             header = FALSE, col.names = c("feature.chr","feature.start","feature.end","signal","bin.chr","bin.start","bin.end")) %>%
mutate(bin = paste0(bin.chr,":",bin.start,"-",bin.end),
       feature_coords = paste0(feature.chr,":",feature.start,"-",feature.end),
       feature = "DNase") %>% 
select(feature, bin, feature_coords, signal) %>% distinct()

head(allregion_dnase,n=1)
dim(allregion_dnase)

allregion_atac <- read.delim("/data/gersbachlab/lrb53/encodeCrisprWG/finalversion/figure1/k562meta/outs/comparesignal/binnedgenome.int.k562.elementReference.int.k562.atac.foldchange.txt",
                             header = FALSE, col.names = c("feature.chr","feature.start","feature.end","signal","bin.chr","bin.start","bin.end")) %>%
mutate(bin = paste0(bin.chr,":",bin.start,"-",bin.end),
       feature_coords = paste0(feature.chr,":",feature.start,"-",feature.end),
       feature = "ATAC") %>% 
select(feature, bin, feature_coords, signal) %>% distinct()

head(allregion_atac,n=1)
dim(allregion_atac)

allregion_h3k4me1 <- read.delim("/data/gersbachlab/lrb53/encodeCrisprWG/finalversion/figure1/k562meta/outs/comparesignal/binnedgenome.int.k562.elementReference.int.k562.h3k4me1.foldchange.txt",
                             header = FALSE, col.names = c("feature.chr","feature.start","feature.end","signal","bin.chr","bin.start","bin.end")) %>%
mutate(bin = paste0(bin.chr,":",bin.start,"-",bin.end),
       feature_coords = paste0(feature.chr,":",feature.start,"-",feature.end),
       feature = "H3K4me1") %>% 
select(feature, bin, feature_coords, signal) %>% distinct()

head(allregion_atac,n=1)
dim(allregion_atac)

allregion_h3k4me3 <- read.delim("/data/gersbachlab/lrb53/encodeCrisprWG/finalversion/figure1/k562meta/outs/comparesignal/binnedgenome.int.k562.elementReference.int.k562.h3k4me3.foldchange.txt",
                             header = FALSE, col.names = c("feature.chr","feature.start","feature.end","signal","bin.chr","bin.start","bin.end")) %>%
mutate(bin = paste0(bin.chr,":",bin.start,"-",bin.end),
       feature_coords = paste0(feature.chr,":",feature.start,"-",feature.end),
       feature = "H3K4me3") %>% 
select(feature, bin, feature_coords, signal) %>% distinct()

head(allregion_atac,n=1)
dim(allregion_atac)

allregion_h3k9me3 <- read.delim("/data/gersbachlab/lrb53/encodeCrisprWG/finalversion/figure1/k562meta/outs/comparesignal/binnedgenome.int.k562.elementReference.int.k562.h3k9me3.foldchange.txt",
                             header = FALSE, col.names = c("feature.chr","feature.start","feature.end","signal","bin.chr","bin.start","bin.end")) %>%
mutate(bin = paste0(bin.chr,":",bin.start,"-",bin.end),
       feature_coords = paste0(feature.chr,":",feature.start,"-",feature.end),
       feature = "H3K9me3") %>% 
select(feature, bin, feature_coords, signal) %>% distinct()

head(allregion_atac,n=1)
dim(allregion_atac)

allregion_h3k27me3 <- read.delim("/data/gersbachlab/lrb53/encodeCrisprWG/finalversion/figure1/k562meta/outs/comparesignal/binnedgenome.int.k562.elementReference.int.k562.h3k27me3.foldchange.txt",
                             header = FALSE, col.names = c("feature.chr","feature.start","feature.end","signal","bin.chr","bin.start","bin.end")) %>%
mutate(bin = paste0(bin.chr,":",bin.start,"-",bin.end),
       feature_coords = paste0(feature.chr,":",feature.start,"-",feature.end),
       feature = "H3K27me3") %>% 
select(feature, bin, feature_coords, signal) %>% distinct()

head(allregion_atac,n=1)
dim(allregion_atac)

allregion_ep300 <- read.delim("/data/gersbachlab/lrb53/encodeCrisprWG/finalversion/figure1/k562meta/outs/comparesignal/binnedgenome.int.k562.elementReference.int.k562.ep300.foldchange.txt",
                             header = FALSE, col.names = c("feature.chr","feature.start","feature.end","signal","bin.chr","bin.start","bin.end")) %>%
mutate(bin = paste0(bin.chr,":",bin.start,"-",bin.end),
       feature_coords = paste0(feature.chr,":",feature.start,"-",feature.end),
       feature = "EP300") %>% 
select(feature, bin, feature_coords, signal) %>% distinct()

head(allregion_atac,n=1)
dim(allregion_atac)

allregion_ctcf <- read.delim("/data/gersbachlab/lrb53/encodeCrisprWG/finalversion/figure1/k562meta/outs/comparesignal/binnedgenome.int.k562.elementReference.int.k562.ctcf.foldchange.txt",
                             header = FALSE, col.names = c("feature.chr","feature.start","feature.end","signal","bin.chr","bin.start","bin.end")) %>%
mutate(bin = paste0(bin.chr,":",bin.start,"-",bin.end),
       feature_coords = paste0(feature.chr,":",feature.start,"-",feature.end),
       feature = "CTCF") %>% 
select(feature, bin, feature_coords, signal) %>% distinct()

head(allregion_atac,n=1)
dim(allregion_atac)

allregion_polr2a <- read.delim("/data/gersbachlab/lrb53/encodeCrisprWG/finalversion/figure1/k562meta/outs/comparesignal/binnedgenome.int.k562.elementReference.int.k562.polr2a.foldchange.txt",
                             header = FALSE, col.names = c("feature.chr","feature.start","feature.end","signal","bin.chr","bin.start","bin.end")) %>%
mutate(bin = paste0(bin.chr,":",bin.start,"-",bin.end),
       feature_coords = paste0(feature.chr,":",feature.start,"-",feature.end),
       feature = "POLR2A") %>% 
select(feature, bin, feature_coords, signal) %>% distinct()

head(allregion_atac,n=1)
dim(allregion_atac)

##### Significant regions

sigregion_h3k27ac <- read.delim("/data/gersbachlab/lrb53/encodeCrisprWG/finalversion/figure1/k562meta/outs/comparesignal/binnedgenome.int.k562.elementQuant.int.k562.h3k27ac.foldchange.txt",
                             header = FALSE, col.names = c("feature.chr","feature.start","feature.end","signal","bin.chr","bin.start","bin.end")) %>%
mutate(bin = paste0(bin.chr,":",bin.start,"-",bin.end)) %>%
select(bin) %>% distinct() %>% pull()

length(sigregion_h3k27ac)

sigregion_dnase <- read.delim("/data/gersbachlab/lrb53/encodeCrisprWG/finalversion/figure1/k562meta/outs/comparesignal/binnedgenome.int.k562.elementQuant.int.k562.dnase.foldchange.txt",
                             header = FALSE, col.names = c("feature.chr","feature.start","feature.end","signal","bin.chr","bin.start","bin.end")) %>%
mutate(bin = paste0(bin.chr,":",bin.start,"-",bin.end)) %>%
select(bin) %>% distinct() %>% pull()

length(sigregion_dnase)

sigregion_atac <- read.delim("/data/gersbachlab/lrb53/encodeCrisprWG/finalversion/figure1/k562meta/outs/comparesignal/binnedgenome.int.k562.elementQuant.int.k562.atac.foldchange.txt",
                             header = FALSE, col.names = c("feature.chr","feature.start","feature.end","signal","bin.chr","bin.start","bin.end")) %>%
mutate(bin = paste0(bin.chr,":",bin.start,"-",bin.end)) %>%
select(bin) %>% distinct() %>% pull()

length(sigregion_atac)

sigregion_h3k4me1 <- read.delim("/data/gersbachlab/lrb53/encodeCrisprWG/finalversion/figure1/k562meta/outs/comparesignal/binnedgenome.int.k562.elementQuant.int.k562.h3k4me1.foldchange.txt",
                             header = FALSE, col.names = c("feature.chr","feature.start","feature.end","signal","bin.chr","bin.start","bin.end")) %>%
mutate(bin = paste0(bin.chr,":",bin.start,"-",bin.end)) %>%
select(bin) %>% distinct() %>% pull()

length(sigregion_atac)

sigregion_h3k4me3 <- read.delim("/data/gersbachlab/lrb53/encodeCrisprWG/finalversion/figure1/k562meta/outs/comparesignal/binnedgenome.int.k562.elementQuant.int.k562.h3k4me3.foldchange.txt",
                             header = FALSE, col.names = c("feature.chr","feature.start","feature.end","signal","bin.chr","bin.start","bin.end")) %>%
mutate(bin = paste0(bin.chr,":",bin.start,"-",bin.end)) %>%
select(bin) %>% distinct() %>% pull()

length(sigregion_atac)

sigregion_h3k9me3 <- read.delim("/data/gersbachlab/lrb53/encodeCrisprWG/finalversion/figure1/k562meta/outs/comparesignal/binnedgenome.int.k562.elementQuant.int.k562.h3k9me3.foldchange.txt",
                             header = FALSE, col.names = c("feature.chr","feature.start","feature.end","signal","bin.chr","bin.start","bin.end")) %>%
mutate(bin = paste0(bin.chr,":",bin.start,"-",bin.end)) %>%
select(bin) %>% distinct() %>% pull()

length(sigregion_atac)

sigregion_h3k27me3 <- read.delim("/data/gersbachlab/lrb53/encodeCrisprWG/finalversion/figure1/k562meta/outs/comparesignal/binnedgenome.int.k562.elementQuant.int.k562.h3k27me3.foldchange.txt",
                             header = FALSE, col.names = c("feature.chr","feature.start","feature.end","signal","bin.chr","bin.start","bin.end")) %>%
mutate(bin = paste0(bin.chr,":",bin.start,"-",bin.end)) %>%
select(bin) %>% distinct() %>% pull()

length(sigregion_atac)

sigregion_ep300 <- read.delim("/data/gersbachlab/lrb53/encodeCrisprWG/finalversion/figure1/k562meta/outs/comparesignal/binnedgenome.int.k562.elementQuant.int.k562.ep300.foldchange.txt",
                             header = FALSE, col.names = c("feature.chr","feature.start","feature.end","signal","bin.chr","bin.start","bin.end")) %>%
mutate(bin = paste0(bin.chr,":",bin.start,"-",bin.end)) %>%
select(bin) %>% distinct() %>% pull()

length(sigregion_atac)

sigregion_ctcf <- read.delim("/data/gersbachlab/lrb53/encodeCrisprWG/finalversion/figure1/k562meta/outs/comparesignal/binnedgenome.int.k562.elementQuant.int.k562.ctcf.foldchange.txt",
                             header = FALSE, col.names = c("feature.chr","feature.start","feature.end","signal","bin.chr","bin.start","bin.end")) %>%
mutate(bin = paste0(bin.chr,":",bin.start,"-",bin.end)) %>%
select(bin) %>% distinct() %>% pull()

length(sigregion_atac)

sigregion_polr2a <- read.delim("/data/gersbachlab/lrb53/encodeCrisprWG/finalversion/figure1/k562meta/outs/comparesignal/binnedgenome.int.k562.elementQuant.int.k562.polr2a.foldchange.txt",
                             header = FALSE, col.names = c("feature.chr","feature.start","feature.end","signal","bin.chr","bin.start","bin.end")) %>%
mutate(bin = paste0(bin.chr,":",bin.start,"-",bin.end)) %>%
select(bin) %>% distinct() %>% pull()

length(sigregion_atac)

plotdf <-
rbind(
    allregion_h3k4me1 %>%
    mutate(result = ifelse(bin %in% sigregion_h3k4me1, "Significant","Not significant")) %>%
    separate(feature_coords, into = c("chr","start","end")),
    allregion_h3k4me3 %>%
    mutate(result = ifelse(bin %in% sigregion_h3k4me1, "Significant","Not significant")) %>%
    separate(feature_coords, into = c("chr","start","end")),
    allregion_h3k9me3 %>%
    mutate(result = ifelse(bin %in% sigregion_h3k9me3, "Significant","Not significant")) %>%
    separate(feature_coords, into = c("chr","start","end")),
    allregion_h3k27ac %>%
    mutate(result = ifelse(bin %in% sigregion_h3k27ac, "Significant","Not significant")) %>%
    separate(feature_coords, into = c("chr","start","end")),
    allregion_h3k27me3 %>%
    mutate(result = ifelse(bin %in% sigregion_h3k27me3, "Significant","Not significant")) %>%
    separate(feature_coords, into = c("chr","start","end"))
    ) %>%
mutate(region_size = as.numeric(end) - as.numeric(start),
       signal_scaled = signal/region_size) %>%
bind_rows(
    rbind(
        allregion_ep300 %>%
        mutate(result = ifelse(bin %in% sigregion_ep300, "Significant","Not significant")) %>%
        separate(feature_coords, into = c("chr","start","end")),
        allregion_ctcf %>%
        mutate(result = ifelse(bin %in% sigregion_ctcf, "Significant","Not significant")) %>%
        separate(feature_coords, into = c("chr","start","end")),
        allregion_polr2a %>%
        mutate(result = ifelse(bin %in% sigregion_polr2a, "Significant","Not significant")) %>%
        separate(feature_coords, into = c("chr","start","end")),
        allregion_atac %>%
        mutate(result = ifelse(bin %in% sigregion_atac, "Significant","Not significant")) %>%
        separate(feature_coords, into = c("chr","start","end")),
        allregion_dnase %>%
        mutate(result = ifelse(bin %in% sigregion_dnase, "Significant","Not significant")) %>%
        separate(feature_coords, into = c("chr","start","end"))
    ) %>%mutate(region_size = as.numeric(end) - as.numeric(start), signal_scaled = signal/region_size)
    )

plotdf$feature = factor(plotdf$feature, levels = c("H3K27me3",
                                                   "H3K9me3",
                                                   "POLR2A",
                                                   "EP300",
                                                   "CTCF",
                                                   "H3K27ac",
                                                   "H3K4me3",
                                                   "H3K4me1",
                                                   "DNase",
                                                   "ATAC"))

rm(list=ls(pattern="allregion_"))
rm(list=ls(pattern="sigregion_"))

str(as.list(.GlobalEnv))

head(plotdf,n=1)

group.colors <- c(`Not significant` = "grey", Significant = "green")

pdf(file="/data/gersbachlab/lrb53/encodeCrisprWG/finalversion/figure1/k562meta/outs/plot_pdfs/signalcomparison.log10transform.scalebyregionsize.k562.pdf")

plotdf %>%
select(feature, signal_scaled, result) %>%
ggplot(aes(y = feature, x = log10(signal_scaled+0.01), fill = result)) +
geom_boxplot() +
ylab("Feature") +
xlab("Log10(Scaled signal) (fold change over background)") +
ggtitle("Feature signal in perturbed genomic regions") +
scale_fill_manual(values = group.colors) +
theme_classic() +
theme(plot.title = element_text(hjust = 0.5),
      legend.position = "bottom"
     )

dev.off()
