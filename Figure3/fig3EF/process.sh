#!/bin/sh
####### Author: Jin Woo Oh  ############
###### For CRISPR WG paper  ############
###### Last updated: 04/16/2023#########
########################################

# sampling for each of the CRISPR screen dataset 

# Bassik GATA1 I
python sampling.py T7_I_Bassik_guide_quantifications.bed  T21_I_R1_Bassik_guide_quantifications.bed T7_I_Bassik_guide_quantifications.bed  T21_I_R2_Bassik_guide_quantifications.bed 100 10 I_Bassik_GATA1

# Bassik GATA1 K
python3 script/sampling.py T7_K_Bassik_guide_quantifications.bed  T21_K_R1_Bassik_guide_quantifications.bed T7_K_Bassik_guide_quantifications.bed  T21_K_R2_Bassik_guide_quantifications.bed 100 10 K_Bassik_GATA1

# Sabeti GATA1  
python script/sampling.py GATA1_HS_exp_r1.tsv GATA1_LS_exp_r1.tsv GATA1_HS_exp_r2.tsv GATA1_LS_exp_r2.tsv 100 10 HCRFF_GATA1

# engreitz GATA1 FF 
python3 script/sampling.py GATA1-1-A-PCR_summed.hg38.tsv   GATA1-1-F-PCR_summed.hg38.tsv  GATA1-3-A-PCR_summed.hg38.tsv   GATA1-3-F-PCR_summed.hg38.tsv 100 10 FF_GATA1


# engreitz GATA1 growth
python3 script/sampling.py doughty_GATA1_growth.t0_rep1.encode_format.hg38.tsv  doughty_GATA1_growth.tend_rep1.encode_format.hg38.tsv doughty_GATA1_growth.t0_rep2.encode_format.hg38.tsv  doughty_GATA1_growth.tend_rep2.encode_format.hg38.tsv  100 10 Engreitz_GATA1_growth

# Sabeti MYC
python script/sampling.py MYC_HS_exp_r1.tsv MYC_LS_exp_r1.tsv MYC_HS_exp_r2.tsv MYC_LS_exp_r2.tsv 100 10 HCRFF_MYC


# engreitz MYC growth 
python3 script/sampling.py doughty_MYC_growth.t0_rep1.encode_format.hg38.tsv  doughty_MYC_growth.tend_rep1.encode_format.hg38.tsv doughty_MYC_growth.t0_rep2.encode_format.hg38.tsv  doughty_MYC_growth.tend_rep2.encode_format.hg38.tsv  100 10 Engreitz_MYC_growth
