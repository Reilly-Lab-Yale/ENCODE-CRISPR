###############################################################################
# Josh Tycko
# 06/05/2023
###############################################################################

'''
Integrate tiling screen data sets

Input: Guide enrichments '_rhos.csv' from casTLE for all perturbations & annotations

Output: Tilng screen plots for encode guideline revision

'''

# Import neccessary modules
import sys
import argparse
import os
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as plticker
import scipy.stats as stats
from scipy.spatial import distance
from scipy.cluster import hierarchy


# mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=["r", "b", "k"]) 
mpl.rcParams['font.sans-serif'] = "Arial"
mpl.rcParams['font.family'] = "sans-serif"
mpl.rcParams['pdf.fonttype'] = 42

###############################################################################    
# Version number

current_version = '0.2'
###############################################################################
# Parses input using argparse module

# Initiates argument parser
parser = argparse.ArgumentParser(description='Plot guide enrichments for tiling screen')

# Non-optional arguments:
parser.add_argument('lib_file', help='File with library annotations',type=str, default = 'Design/TilingLibrary.xlsx')

parser.add_argument('KR1_file', help='File with KR1 enrichments',
                    type=str, default = 'Data/T14vsPlasmid_K_R1_Tiling_rhos.csv')

parser.add_argument('KR2_file', help='File with KR2 enrichments',
                    type=str, default = 'Data/T14vsPlasmid_K_R2_Tiling_rhos.csv')

parser.add_argument('IR1_file', help='File with IR1 enrichments',
                    type=str, default = 'Data/T14vsPlasmid_I_R1_Tiling_rhos.csv')

parser.add_argument('IR2_file', help='File with IR2 enrichments',
                    type=str, default = 'Data/T14vsPlasmid_I_R2_Tiling_rhos.csv')

parser.add_argument('AR1_file', help='File with AR1 enrichments',
                    type=str, default = 'Data/T14vsPlasmid_A_R1_Tiling_rhos.csv')

parser.add_argument('AR2_file', help='File with AR2 enrichments',
                    type=str, default = 'Data/T14vsPlasmid_A_R2_Tiling_rhos.csv')

parser.add_argument('DR1_file', help='File with DR1 enrichments',
                    type=str, default = 'Data/T14vsPlasmid_D_R1_Tiling_rhos.csv')

parser.add_argument('DR2_file', help='File with DR2 enrichments',
                    type=str, default = 'Data/T14vsPlasmid_D_R2_Tiling_rhos.csv')


# Optional arguments:
parser.add_argument('-p', '--plotType', dest='plotType',
                    help='The type of plot to generate', choices = ['gene', 'specificity'],
                    type=str, default='specificity')

# Saves input to args object
args = parser.parse_args()

###############################################################################
# # Makes checks for failed arguments

# # Determines output file
# file_out = os.path.join('Plots',args.name)

# try:
#     with open(file_out + '_GuideEnrichments.pdf','w') as file_open:
#         pass
# except:
#     sys.exit('Cannot write to output file: \n' + file_out + '\n')

###############################################################################
# Read input files and create dataframe

libdf = pd.read_excel(args.lib_file, engine='openpyxl')

KR1df = pd.read_csv(args.KR1_file, names = ['label','KR1'], header = None)
KR2df = pd.read_csv(args.KR2_file, names = ['label','KR2'], header = None)

AR1df = pd.read_csv(args.AR1_file, names = ['label','AR1'], header = None)
AR2df = pd.read_csv(args.AR2_file, names = ['label','AR2'], header = None)

IR1df = pd.read_csv(args.IR1_file, names = ['label','IR1'], header = None)
IR2df = pd.read_csv(args.IR2_file, names = ['label','IR2'], header = None)

DR1df = pd.read_csv(args.DR1_file, names = ['label','DR1'], header = None)
DR2df = pd.read_csv(args.DR2_file, names = ['label','DR2'], header = None)

# df = pd.merge(libdf, BR1df, on = 'label', how = 'outer')
df = KR1df
for dfx in [KR2df, AR1df , AR2df, IR1df , IR2df, DR1df, DR2df]:
	df = pd.merge(df, dfx, on = 'label', how = 'outer')

# df['Avg'] = df[['BR1','BR2']].mean(axis = 1)
# df['Standard Error'] = df[['BR1','BR2']].sem(axis = 1)

df['Gene'],df['Location'],df['blank'],df['Guide ID'] = df['label'].str.split('_').str
df = df.drop(['blank'], axis = 1)
df['Chromosome'] = df['Location'].str.split(':').str[0]
df['Strand'] = df['Location'].str.split(':').str[-1]
df['Guide start'], df['Guide end'] = df['Location'].str.split(':').str[1].str.split('-').str

def find_cut_site(row):
	if row['Strand'] == '-':
		return int(row['Guide start']) + 3
	elif row['Strand'] == '+':
		return int(row['Guide end']) - 3

df['Cut site'] = df.apply(find_cut_site, axis = 1)
df['Guide start'] = pd.to_numeric(df['Guide start'] , errors = 'coerce')
df = pd.merge(df, libdf, on = 'Guide start', how = 'inner')


# elevdf = pd.read_table(args.Elev_file, sep = " ", names = ['Guide_PAM', 'Elevation-aggregate'])
# elevdf['guide'] = elevdf['Guide_PAM'].str[:-4]
# df = pd.merge(df, elevdf, on = 'guide', how = 'outer')

###############################################################################
if args.plotType == 'gene':
	gene = 'GATA1'
	GATAdf = df.loc[df['Gene'] == gene]

	eGATA1_df = GATAdf.loc[GATAdf['Guide start'].isin(range(48641136,48641797))]
	eHDAC6_df = GATAdf.loc[GATAdf['Guide start'].isin(range(48658755,48658755+700))]
	TSS_df = GATAdf.loc[GATAdf['Guide start'].isin(range(48644981-200,48644981+200))]
	exon_df = GATAdf.loc[GATAdf['Growth gene exon'] > 0]
	CFD_df = GATAdf.loc[(GATAdf['cutting specificity score'] < 0.2) & (GATAdf[['DR1','DR2','KR1','KR2', 'IR1','IR2', 'AR1','AR2']].min(axis = 1) < -3) & (GATAdf['Growth gene exon'] == 0) & (~GATAdf['Guide start'].isin(range(48641136-1000,48658755+1700)))] # Low CFD guides with effect, outside other annotations

	dfs = ['eGATA1' , 'eHDAC6', 'TSS', 'exon', 'Low_CFD']
	idx = 0
	for dfx in [eGATA1_df, eHDAC6_df, TSS_df, exon_df, CFD_df]:
		D = distance.pdist(dfx[['DR1','DR2','KR1','KR2', 'IR1','IR2', 'AR1','AR2']])
		# Dt = distance.pdist(df.T)
		row_linkage = hierarchy.linkage(D, method='average')
		# col_linkage = hierarchy.linkage(Dt, method='average')
		optimal_Row = optimal_leaf_ordering(row_linkage, D)
		# optimal_Col = optimal_leaf_ordering(col_linkage, Dt)
		g = sns.clustermap(data = dfx[['DR1','DR2','KR1','KR2', 'IR1','IR2', 'AR1','AR2']], row_linkage = row_linkage, vmin = -4, vmax = 1, center =0, yticklabels = True, col_cluster = False, cmap = "RdBu_r")
		# plt.show()
		plt.savefig('../Figures/Cluster/' + gene + '_cluster_' + dfs[idx] + '.pdf', transparent=True)
		plt.close()
		idx += 1



# Koff_df = df.loc[(df['cutting specificity score'] < 0.2) & (df[['KR1','KR2']].min(axis = 1) < -3) & (df['Growth gene exon'] == 0) & (~df['Guide start'].isin(range(48641136-2000,48658755+2700)))] # Low CFD guides with effect, outside other annotations
# print('KO low CFD guides', 'aggn'.join(Koff_df['guide']))
# Ioff_df = df.loc[(df['cutting specificity score'] < 0.2) & (df[['IR1','IR2']].min(axis = 1) < -3) & (df['Growth gene exon'] == 0) & (~df['Guide start'].isin(range(48641136-2000,48658755+2700)))] # Low CFD guides with effect, outside other annotations
# print('I low CFD guides', 'aggn'.join(Ioff_df['guide']))
# Aoff_df = df.loc[(df['cutting specificity score'] < 0.2) & (df[['AR1','AR2']].min(axis = 1) < -3) & (df['Growth gene exon'] == 0) & (~df['Guide start'].isin(range(48641136-2000,48658755+2700)))] # Low CFD guides with effect, outside other annotations
# print('A low CFD guides', 'aggn'.join(Aoff_df['guide']))
# CtrlOff_df = df.loc[(df['cutting specificity score'] < 0.2) & (abs(df[['DR1','DR2','KR1','KR2', 'IR1','IR2', 'AR1','AR2']]).max(axis = 1) < .333) & (df['Growth gene exon'] == 0) & (~df['Guide start'].isin(range(48641136-2000,48658755+2700)))]
# print('Control low CFD guides with no effect', 'aggn'.join(CtrlOff_df['guide']))


df['K_Avg'] = df[['KR1','KR2']].mean(axis = 1)
df['A_Avg'] = df[['AR1','AR2']].mean(axis = 1)
df['I_Avg'] = df[['IR1','IR2']].mean(axis = 1)
df['D_Avg'] = df[['DR1','DR2']].mean(axis = 1)
df['offtargets sum +1 log']=np.log10(df['offtargets sum'] +1)


eGATA1_df = df.loc[df['Guide start'].isin(range(48641136,48641797))]
eHDAC6_df = df.loc[df['Guide start'].isin(range(48658755,48658755+700))]
GATA1_TSS_df = df.loc[df['Guide start'].isin(range(48644981-250,48644981+250))]
MYB_TSS_df = df.loc[df['Guide start'].isin(range(135502452-250,135502452+250))]
ZMYND8_TSS_df = df.loc[df['Guide start'].isin(range(45985633-250,45985633+250))]

enhTSS_df = pd.concat([eGATA1_df, eHDAC6_df, GATA1_TSS_df, MYB_TSS_df, ZMYND8_TSS_df])

exon_df = df.loc[(df['Growth gene exon'] > 0)]
g = sns.jointplot(x ='cutting specificity score', y = 'K_Avg', data = exon_df, alpha = 0.3, linewidth = 0, space = 0, ylim = (-6.2, 2.2), xlim = (-.05, 1.05))
g.set_axis_labels('GuideScan CFD Specificity Score', 'Guide enrichment')
plt.savefig('Figures/Off-targets/T14vsPlasmid_K_GrowthExons_CFD_Model.pdf', transparent = True)
# plt.show()
plt.close()

# g = sns.jointplot(x ='cutting specificity score', y = 'I_Avg', data = enhTSS_df, alpha = 0.3, linewidth = 0, space = 0, ylim = (-6.2, 2.2), xlim = (-.05, 1.05))
# g.set_axis_labels('GuideScan CFD Specificity Score', 'Guide enrichment')
# plt.savefig('../Figures/Off-targets/T14vsPlasmid_I_GrowthEnhTSS_CFD_Model.pdf', transparent = True)
# # plt.show()
# plt.close()

# g = sns.jointplot(x ='cutting specificity score', y = 'A_Avg', data = enhTSS_df, alpha = 0.3, linewidth = 0, space = 0, ylim = (-6.2, 2.2), xlim = (-.05, 1.05))
# g.set_axis_labels('GuideScan CFD Specificity Score', 'Guide enrichment')
# plt.savefig('../Figures/Off-targets/T14vsPlasmid_A_GrowthEnhTSS_CFD_Model.pdf', transparent = True)
# # plt.show()
# plt.close()

#filter for likley Off targets
def in_region(chrom, start, end):
    return ((df['Chromosome_x'] == chrom) &
            (start <= df['Guide start']) & (df['Guide start'] < end))

df['likely_on_target'] = (
    in_region('chrX', 48641136 - 1000, 48659966 + 1000) |
    in_region('chr20', 45837859 - 1000, 45985567 + 1000) |
    in_region('chr6', 135502453 - 1000, 135540311 + 1000)
)
subdf = df[~df['likely_on_target']]

# focus on GATA1 tiling
# subdf = subdf.loc[subdf['Chromosome_x'] == 'chrX']
# exon_df = exon_df.loc[exon_df['Chromosome_x'] == 'chrX']
# enhTSS_df = enhTSS_df.loc[enhTSS_df['Chromosome_x'] == 'chrX']
for perturb in ['K','I','A','D']:
	g = sns.jointplot(x ='cutting specificity score', y = perturb+'_Avg', data = subdf, alpha = 0.3, linewidth = 0, space = 0, ylim = (-6.2, 2.2), xlim = (-.05, 1.05))
	# if perturb == 'K':
	# 	g.x = exon_df['cutting specificity score']
	# 	g.y = exon_df['K_Avg']
	# 	g.plot_joint(plt.scatter, c='grey')
	# if perturb == 'I':
	# 	g.x = enhTSS_df['cutting specificity score']
	# 	g.y = enhTSS_df['I_Avg']
	# 	g.plot_joint(plt.scatter, c='grey')
	# # if perturb == 'A':
	# # 	g.x = GATA1_TSS_df['cutting specificity score']
	# # 	g.y = GATA1_TSS_df['A_Avg']
	# # 	g.plot_joint(plt.scatter, c='grey')
	g.set_axis_labels('GuideScan specificity score', 'Guide enrichment')
	plt.savefig('Figures/Off-targets/T14vsPlasmid_'+perturb+'_offTarget_GuideScans.pdf', transparent = True, rasterize = True)
	# plt.show()
	plt.close()
	top_left = sum((subdf['cutting specificity score'] < 0.2) & (subdf[perturb+'_Avg'] > -2.0))
	bottom_left = sum((subdf['cutting specificity score'] < 0.2) & (subdf[perturb+'_Avg'] <= -2.0))
	top_right = sum((subdf['cutting specificity score'] >= 0.2) & (subdf[perturb+'_Avg'] > -2.0))
	bottom_right = sum((subdf['cutting specificity score'] >= 0.2) & (subdf[perturb+'_Avg'] <= -2.0))
	# sum([top_right, top_left, bottom_right, bottom_left])
	print('Quadrants in GuideScan chart for CRISPR_' + perturb, [top_left, top_right, bottom_left, bottom_right], stats.fisher_exact([[top_left, top_right], [bottom_left, bottom_right]]))
	# OT count plots
	g = sns.jointplot(x ='offtargets sum +1 log', y = perturb + '_Avg', data = df, alpha = 0.3, linewidth = 0, space = 0, ylim = (-6.2, 2.2))
	g.set_axis_labels('log10(Count of off-targets +1)', 'Guide enrichment')
	plt.savefig('Figures/Off-targets/T14vsPlasmid_'+perturb+'_LogCountOT.pdf', transparent = True, rasterize = True)
	###  redo this with subsetted for outside the ontarget region
	g = sns.jointplot(x ='offtargets sum +1 log', y = perturb + '_Avg', data = subdf, alpha = 0.3, linewidth = 0, space = 0, ylim = (-6.2, 2.2))
	g.set_axis_labels('log10(Count of off-targets +1)', 'Guide enrichment')
	plt.savefig('Figures/Off-targets/T14vsPlasmid_'+perturb+'_LogCountOT_likelyOFF.pdf', transparent = True, rasterize = True)

# GuideScan score cutoff optimization
for perturb in ['K','I','A','D']:
	totalBad_likelyOFF = sum(subdf[perturb+'_Avg']<-2)
	x  = []
	y1 = []
	y2 = []
	for c in np.arange(0,1.05,0.05):
		x.append(c)
		y1.append(1-(subdf.loc[(subdf[perturb+'_Avg']<-2) & (subdf['cutting specificity score']>=c)]).shape[0]/float(totalBad_likelyOFF))
		y2.append(sum(df['cutting specificity score']>=c)/float(df.shape[0]))
	
	fig, ax1 = plt.subplots()
	color = 'tab:red'
	ax1.set_xlabel('GuideScan specificity score cutoff')
	ax1.set_ylabel('Fraction confounded sgRNAs removed', color=color)
	ax1.set_ylim([0, 1.1])
	ax1.plot(x, y1, color=color)
	ax1.tick_params(axis='y', labelcolor=color)
	#
	ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
	color = 'tab:blue'
	ax2.set_ylabel('Fraction library remaining', color=color)  # we already handled the x-label with ax1
	ax2.set_ylim([0, 1.1])
	ax2.plot(x, y2, color=color)
	ax2.tick_params(axis='y', labelcolor=color)
	#
	fig.tight_layout()  # otherwise the right y-label is slightly clipped
	plt.savefig('Figures/Off-targets/T14vsPlasmid_'+perturb+'_GScutoffOpt_likelyOFF.pdf', transparent = True, rasterize = True)
	plt.close()


# Cas9 vs CRISPRi at promoters and enhancers
def labelTarget(row):
     if row['Guide start'] in range(48641136,48641797):
          target = 'eGATA1'
     elif row['Guide start'] in range(48658755,48658755+700):
          target = 'eHDAC6'
     elif row['Guide start'] in range(48644981-500,48644981+500):
          target = 'GATA1 TSS'
     elif row['Guide start'] in range(135502452-500,135502452+500):
          target = 'MYB TSS'
     elif row['Guide start'] in range(45985633-500,45985633+500):
          target = 'ZMYND8 TSS'
     elif row['Gene']=='GATA1' and row['Growth gene exon']>0:
          target = 'GATA1 exon'
     else:
          target = 'other sgRNA'
     return target

df['sgRNA target'] = df.apply(labelTarget, axis = 1)

def labelExon(row):
    if row['Gene']=='GATA1' and row['Growth gene exon']>0:
        exon = True
    else:
     	exon = False
    return exon

df['GATA1 exon'] = df.apply(labelExon, axis = 1)


tempdf = df.loc[(df['cutting specificity score']>0.2) & (df['sgRNA target'].isin(['GATA1 exon','eGATA1','eHDAC6','GATA1 TSS']))]
tempdf['I_Avg'] = -1*tempdf['I_Avg']
tempdf['A_Avg'] = -1*tempdf['A_Avg']
tempdf['K_Avg'] = -1*tempdf['K_Avg']
tempdf['D_Avg'] = -1*tempdf['D_Avg']
g = sns.scatterplot(x = 'I_Avg', y = 'K_Avg' , data = tempdf, hue='sgRNA target', style='GATA1 exon', palette=['pink','k','r','darkred'], s=100, legend = True)
g.set_xlabel('CRISPRi log2FC')
g.set_ylabel('Cas9 log2FC')
plt.savefig('Figures/Modality/Cas9vsI_GrowthEnhTSS_CFDFilt.pdf', transparent = True)
plt.close() 
#CRISPRi vs CRISPRa in same way
g = sns.scatterplot(x = 'I_Avg', y = 'A_Avg' , data = tempdf, hue='sgRNA target', style='GATA1 exon', palette=['pink','k','r','darkred'], s=100, legend = True)
g.set_xlabel('CRISPRi log2FC') 
g.set_ylabel('CRISPRa log2FC')
plt.savefig('Figures/Modality/CRISPRavsI_GrowthEnhTSS_CFDFilt.pdf', transparent = True)
plt.close() 
#CRISPRd vs CRISPRa in same way
g = sns.scatterplot(x = 'D_Avg', y = 'A_Avg' , data = tempdf, hue='sgRNA target', style='GATA1 exon', palette=['pink','k','r','darkred'], s=100, legend = True)
g.set_xlabel('dCas9 log2FC') 
g.set_ylabel('CRISPRa log2FC')
plt.savefig('Figures/Modality/CRISPRavsdCas9_GrowthEnhTSS_CFDFilt.pdf', transparent = True)
plt.close() 




