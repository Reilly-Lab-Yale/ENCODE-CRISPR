#!/bin/sh
####### Author: Jin Woo Oh  ############
###### For CRISPR WG paper  ############
###### Last updated: 04/16/2023#########
########################################

for c in 20 50 100 200
do
	
	for s in 1 5  10 20 50 100 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900 2000
	do
	
	# bootstrap sample gRNA mapped reads
	python sampling_generate_guidequant.py  GATA_${c}x_HS_ENCODE_guideQuant.tsv ${s} 10 GATA_${c}x_HS_ENCODE_quideQuant_bootstrap_${s}x
	python sampling_generate_guidequant.py  GATA_${c}x_LS_ENCODE_guideQuant.tsv ${s} 10 GATA_${c}x_LS_ENCODE_quideQuant_bootstrap_${s}x


	# compute log2FC and annotate as CRE overlapping or not 
	python3 gRNA_to_log2FC.py  GATA_${c}x_LS_ENCODE_quideQuant_bootstrap_${s}x_sim_${i}_guidequant.txt GATA_${c}x_HS_ENCODE_quideQuant_bootstrap_${s}x_sim_${i}_guidequant.txt  GATA_${c}x_ENCODE_quideQuant_bootstrap_${s}x_sim_${i} 
	python3 summary_to_bedgraph.py  GATA_${c}x_ENCODE_quideQuant_bootstrap_${s}x_sim_${i}.log2FC_summary 5 GATA_${c}x_ENCODE_quideQuant_bootstrap_${s}x_sim_${i}.log2FC_3.bedgraph
	python annotate_by_CASA_peak_overlap.py   GATA_${c}x_ENCODE_quideQuant_bootstrap_${s}x_sim_${i}.log2FC_3.bedgraph  ENCFF413WYU.bed > GATA_${c}x_ENCODE_quideQuant_bootstrap_${s}x_sim_${i}.log2FC_3.bedgraph.peak_overlap


	# compute AUPRC
	python compute_auprc.py GATA_${c}x_ENCODE_quideQuant_boostrap_${s}x_sim_${i}.log2FC_3.bedgraph.peak_overlap 100 ${c} ${s} ${i} cc_${c}_sd_${s}_sid_${i}.auprc

	done
done



# generate table for plot
printf "cellcov\tseqdepth\tsample_id\tauprc\n" > all_cellcov_seqdepth_bootstrapsample.aurpc
for file in *.auprc
do
        cat $file >> all_cellcov_seqdepth_bootstrapsample.aurpc
        printf "\n" >> all_cellcov_seqdepth_bootstrapsample.aurpc
done

