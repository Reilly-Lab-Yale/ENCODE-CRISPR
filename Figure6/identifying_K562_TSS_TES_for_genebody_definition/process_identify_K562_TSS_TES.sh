#!/bin/sh
####### Author: Jin Woo Oh  ############
###### For CRISPR WG paper  ############
###### Last updated: 04/23/2023#########
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
python3 index_to_genebody.py  K562_best_isoform_indices.txt  refGene.txt|sort|uniq > K562_gene_bodies_selected_by_higest_Pol2_TSS_TES.out
# for K562 HCRFF covered genes, the above  unambiguously provided us the highest confidence gene bodies. 


# not used in the crispr paper -- we can also generate the best isoform list using TSS pol2 signal only. 
python index_to_TSS.py isoform_index_with_K562_pol2_TSS.out   refGene.txt |sort|uniq > K562_TSS_selected_by_higest_Pol2_TSS.out
python annotate_with_pol2_signal.py  > K562_TSS_selected_by_higest_Pol2_TSS_sig.out
