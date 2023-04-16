#!/bin/sh
####### Author: Jin Woo Oh  ############
###### For CRISPR WG paper  ############
###### Last updated: 04/16/2023#########
########################################


## this shell script is used to compute strand specific log2FC for Reily/Sabeti lab HCRFF, but the same pipeline is used for other CRISPR dataset.

# FADS1 gene coordinate
gtf="refGene.txt"
dhs="ENCFF274YGF.bed"

# strand specific
for gene in CAPRIN1 CAT     CD164   ERP29   FADS1   FADS2   FADS3   FEN1    GATA1   HBE1    HBG1    HBG2    HBS1L   HDAC6   LMO2    MEF2C   MYB     MYC     NMU     PVT1
do

	# For bioreplicate 1&2, compute log2FC separately for guides targeting '+' and '-' strand
	for i in 1 2 
	do
		head -1 ${gene}_HCRFF_r${i}.log2FC_summary > ${gene}_HCRFF_pstrand_r${i}.log2FC_summary
		cat ${gene}_HCRFF_r${i}.log2FC_summary|grep :+ >> ${gene}_HCRFF_pstrand_r${i}.log2FC_summary

		head -1 ${gene}_HCRFF_r${i}.log2FC_summary > ${gene}_HCRFF_nstrand_r${i}.log2FC_summary
		cat ${gene}_HCRFF_r${i}.log2FC_summary|grep :- >> ${gene}_HCRFF_nstrand_r${i}.log2FC_summary

		# the same summary_to_bedgraph.py used in fig3
		python3 summary_to_bedgraph.py  ${gene}_HCRFF_pstrand_r${i}.log2FC_summary ztransformed_by_all_guides ${gene}_HCRFF_pstrand_r${i}.log2FC_Z.bedgraph
		python3 summary_to_bedgraph.py  ${gene}_HCRFF_nstrand_r${i}.log2FC_summary ztransformed_by_all_guides ${gene}_HCRFF_nstrand_r${i}.log2FC_Z.bedgraph

	done

	# average the two biorep log2FC
	python3 avg_reps.py  ${gene}_HCRFF_nstrand_r1.log2FC_Z.bedgraph ${gene}_HCRFF_nstrand_r2.log2FC_Z.bedgraph > ${gene}_HCRFF_nstrand_ravg.log2FC_Z.bedgraph
	python3 avg_reps.py  ${gene}_HCRFF_pstrand_r1.log2FC_Z.bedgraph ${gene}_HCRFF_pstrand_r2.log2FC_Z.bedgraph > ${gene}_HCRFF_pstrand_ravg.log2FC_Z.bedgraph

	#re-label "+" and "-" as template or coding strand by looking up the gene annotatio .gtf file. 
	python3 template_or_coding.py ${gene}_HCRFF_nstrand_ravg.log2FC_Z.bedgraph ${gene} ${gtf}
	python3 template_or_coding.py ${gene}_HCRFF_pstrand_ravg.log2FC_Z.bedgraph ${gene} ${gtf}


	# Using the the K562 gene body annotation, further partition the guides into gene body, promoter, .., -targeting guides. 
	python3 partition_bedgraph_by_gene_body.py ${gene}_HCRFF_template_ravg.log2FC_Z.bedgraph  ${gtf}  ${gene} K562_best_isoform_indices.txt ${dhs} ${gene}_HCRFF_template_gene_body_ravg.log2FC_Z.bedgraph  ${gene}_HCRFF_template_promoter_ravg.log2FC_Z.bedgraph  ${gene}_HCRFF_template_outside_GB_prom_DHS_ravg.log2FC_Z.bedgraph ${gene}_HCRFF_template_TSS_5kb_flank_ravg.log2FC_Z.bedgraph
	python3 partition_bedgraph_by_gene_body.py ${gene}_HCRFF_coding_ravg.log2FC_Z.bedgraph  ${gtf}  ${gene} K562_best_isoform_indices.txt ${dhs} ${gene}_HCRFF_coding_gene_body_ravg.log2FC_Z.bedgraph  ${gene}_HCRFF_coding_promoter_ravg.log2FC_Z.bedgraph ${gene}_HCRFF_coding_outside_GB_prom_DHS_ravg.log2FC_Z.bedgraph ${gene}_HCRFF_coding_TSS_5kb_flank_ravg.log2FC_Z.bedgraph

done
