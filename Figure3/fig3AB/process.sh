#!/bin/sh

####### Author: Jin Woo Oh  ############
###### For CRISPR WG paper  ############
###### Last updated: 04/16/2023#########
########################################

# compute log2FC
#for gene in CAPRIN1 CAT CD164 ERP29 FADS1 FADS2 FADS3 FEN1 GATA1 LMO2 MEF2C MYC NMU PVT1
for i in 20 50 100 200 
do
echo $gene

# compute log2FC
python3 gRNA_to_log2FC.py  GATA_${i}x_LS_ENCODE_guideQuant.tsv GATA_${i}x_HS_ENCODE_guideQuant.tsv GATA_HCRFF_x${i}
python3 summary_to_bedgraph.py  GATA_HCRFF_x${i}.log2FC_summary normalized_by_mean_count1 GATA_HCRFF_x${i}.log2FC_MN1.bedgraph

# overlap with CASA peaks to annotate by CRE 
python annotate_by_CASA_peak_overlap.py  GATA_HCRFF_x${i}.log2FC_MN1.bedgraph  ENCFF413WYU.bed > GATA_HCRFF_x${i}.log2FC_MN1.bedgraph.peak_overlap

# compute ROC
python auroc_auprc.py  GATA_HCRFF_x${i}.log2FC_MN1.bedgraph.peak_overlap 100 precision_recall_x${i}.out

done
