#!/bin/sh
####### Author: Jin Woo Oh  ############
###### For CRISPR WG paper  ############
###### Last updated: 04/16/2023#########
########################################


# use TSS pol2 signal for picking the k562 isoform
# if multiple isoforms with same pol2, also use TES pol2

# TSS
# extract coordinates for TSS flakning 1kb loci
python extract_TSS_flank_1k_bed.py refGene.txt | awk '$2>0' > all_isoform_TSS_flank_1kb_f.bed
# compute avg POL2A chip-seq signals
/kent/bigWigAverageOverBed  ENCFF914WIS.bigWig  all_isoform_TSS_flank_1kb_f.bed all_isoform_TSS_flank_1kb_Pol2_chip.out 
# extract transcript isoform with top pol2 signal
python extract_genefile_index_for_highest_pol2.py  all_isoform_TSS_flank_1kb_Pol2_chip.out > isoform_index_with_K562_pol2_TSS.out


# TES -- same but for TES 
python extract_TES_flank_1k_bed.py refGene.txt | awk '$2>0' > all_isoform_TES_flank_1kb_f.bed
/kent/bigWigAverageOverBed  ENCFF914WIS.bigWig  all_isoform_TES_flank_1kb_f.bed all_isoform_TES_flank_1kb_Pol2_chip.out
python extract_genefile_index_for_highest_pol2.py  all_isoform_TES_flank_1kb_Pol2_chip.out > isoform_index_with_K562_pol2_TES.out



# combining TSS TES pol2 for best isoforms
python combine_TSS_TES_best.py isoform_index_with_K562_pol2_TSS.out isoform_index_with_K562_pol2_TES.out > K562_best_isoform_indices.txt



# for K562 HCRFF covered genes, this method unambiguously provided us the highest confidence gene bodies. 


# to generate list of K562 genes w/ their coordinates. Genes with multiple gene coordinates with equally high TSS/TES Pol2 enrichment are duplicated. 
python3 index_to_TSSfile.py  K562_best_isoform_indices.txt  /mnt/t/data0/joh27/projects/ENCODE/ref_genome/gene_annotation/refGene.txt|uniq > K562_gene_bodies_selected_by_higest_Pol2_TSS_TES.out
