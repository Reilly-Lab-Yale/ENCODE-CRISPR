###############################################################################
# Josh Tycko
# 11/14/2018
# Updating 2/19/2021

'''
python -i ../scripts/18.11.14_TilingScreenTimepoints.py Design/TilingLibrary.xlsx Data/T14vsPlasmid_I_R1_Tiling_rhos.csv Data/T14vsPlasmid_I_R2_Tiling_rhos.csv Data/T14vsT0_I_R1_Tiling_rhos.csv Data/T14vsT0_I_R2_Tiling_rhos.csv Data/T0vsPlasmid_I_Tiling_rhos.csv Data/T14vsPlasmid_D_R1_Tiling_rhos.csv Data/T14vsPlasmid_D_R2_Tiling_rhos.csv
'''
###############################################################################

'''
Integrate tiling screen data sets from two diff timepoints

Input: Guide enrichments '_rhos.csv' from casTLE for 

Output: Plots

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
# from polo import optimal_leaf_ordering


# mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=["r", "b", "k"]) 
mpl.rcParams['font.sans-serif'] = "Arial"
mpl.rcParams['font.family'] = "sans-serif"
sns.set(style="ticks")
mpl.rcParams['pdf.fonttype'] = 42

###############################################################################    
# Version number

current_version = '0.1'
###############################################################################
# Parses input using argparse module

# Initiates argument parser
parser = argparse.ArgumentParser(description='Plot guide enrichments across a region')

# Non-optional arguments:
parser.add_argument('lib_file', help='File with library annotations',
                    type=str)

parser.add_argument('Pl_R1_file', help='File with Plasmid vs T14 R1 enrichments',
                    type=str)

parser.add_argument('Pl_R2_file', help='File with Plasmid vs T14 R2 enrichments',
                    type=str)

parser.add_argument('T0_R1_file', help='File with T0 vs T14 R1 enrichments',
                    type=str)

parser.add_argument('T0_R2_file', help='File with T0 vs T14 R2 enrichments',
                    type=str)

parser.add_argument('T0_Pl_file', help='File with Plasmid vs T0 enrichments (one rep)',
                    type=str)

parser.add_argument('Pl_R1_d_file', help='File with Plasmid vs T14 R1 enrichments for dCas9',
                    type=str)

parser.add_argument('Pl_R2_d_file', help='File with Plasmid vs T14 R2 enrichments for dCas9',
                    type=str)

# python -i ../scripts/18.11.14_TilingScreenTimepoints.py Design/TilingLibrary.xlsx Data/T14vsPlasmid_I_R1_Tiling_rhos.csv Data/T14vsPlasmid_I_R2_Tiling_rhos.csv Data/T14vsT0_I_R1_Tiling_rhos.csv Data/T14vsT0_I_R2_Tiling_rhos.csv Data/T0vsPlasmid_I_Tiling_rhos.csv Data/T14vsPlasmid_D_R1_Tiling_rhos.csv Data/T14vsPlasmid_D_R2_Tiling_rhos.csv

# Optional arguments:
parser.add_argument('-p', '--plotType', dest='plotType',
                    help='The type of plot to generate', choices = ['gene', 'specificity'],
                    type=str, default='region-category')


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

Pl_R1df = pd.read_csv(args.Pl_R1_file, names = ['label','Pl_R1'], header = None)
Pl_R2df = pd.read_csv(args.Pl_R2_file, names = ['label','Pl_R2'], header = None)

T0_R1df = pd.read_csv(args.T0_R1_file, names = ['label','T0_R1'], header = None)
T0_R2df = pd.read_csv(args.T0_R2_file, names = ['label','T0_R2'], header = None)

T0_Pldf = pd.read_csv(args.T0_Pl_file, names = ['label','T0_Pl'], header = None)

Pl_R1df_d = pd.read_csv(args.Pl_R1_d_file, names = ['label','Pl_R1_dCas9'], header = None)
Pl_R2df_d = pd.read_csv(args.Pl_R2_d_file, names = ['label','Pl_R2_dCas9'], header = None)

# df = pd.merge(libdf, BR1df, on = 'label', how = 'outer')
df = Pl_R1df
for dfx in [Pl_R2df, T0_R1df , T0_R2df, T0_Pldf, Pl_R1df_d, Pl_R2df_d]:
     df = pd.merge(df, dfx, on = 'label', how = 'outer')


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

df['Plasmid vs T21 enrichment'] = df[['Pl_R1','Pl_R2']].mean(axis = 1)
df['Plasmid vs T21 error'] = df[['Pl_R1','Pl_R2']].sem(axis = 1)

df['T07 vs T21 enrichment'] = df[['T0_R1','T0_R2']].mean(axis = 1)
df['T07 vs T21 error'] = df[['T0_R1','T0_R2']].sem(axis = 1)

df['Plasmid vs T07 enrichment'] = df['T0_Pl']

df['Plasmid vs T21 dCas9 enrichment'] = df[['Pl_R1_dCas9','Pl_R2_dCas9']].mean(axis = 1)
df['Plasmid vs T21 dCas9 error'] = df[['Pl_R1_dCas9','Pl_R2_dCas9']].sem(axis = 1)


eGATA1_df = df.loc[df['Guide start'].isin(range(48641136,48641797))]
eHDAC6_df = df.loc[df['Guide start'].isin(range(48658755,48658755+700))]
GATA1_TSS_df = df.loc[df['Guide start'].isin(range(48644981-500,48644981+500))]
MYB_TSS_df = df.loc[df['Guide start'].isin(range(135502452-500,135502452+500))]
ZMYND8_TSS_df = df.loc[df['Guide start'].isin(range(45985633-500,45985633+500))]

enhTSS_df = pd.concat([eGATA1_df, eHDAC6_df, GATA1_TSS_df, MYB_TSS_df, ZMYND8_TSS_df])

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
     else:
          target = 'other sgRNA'
     return target


df['sgRNA target'] = df.apply(labelTarget, axis = 1)

###############################################################################
# Make plots

g = sns.pairplot(df[['Pl_R1', 'Pl_R2', 'T0_R1', 'T0_R2','T0_Pl']].dropna())
plt.savefig('Figures/Timepoints/CompareEnrichments.pdf', transparent=True)
plt.close()


g = sns.jointplot(x = 'Plasmid vs T21 enrichment', y = 'T07 vs T21 enrichment' , data = df, alpha = 0.3, linewidth = 0, space = 0)
plt.savefig('Figures/Timepoints/CompareEnrichmentsCRISPRi_Avg.pdf', transparent=True)
plt.close()

data = df.loc[df['Chromosome_x']=='chrX']
for timepoint in ['T07 vs T21 enrichment','Plasmid vs T07 enrichment']:
     sns.lmplot(y = 'Plasmid vs T21 enrichment', x = timepoint , data = data,  hue = 'sgRNA target', fit_reg=False, palette=['tab:gray','k','b','r','y','m'], legend = True, scatter_kws={'alpha':0.7,'marker':'.'})
     colordict = {'other sgRNA':'tab:gray','GATA1 TSS':'k','eHDAC6':'b','eGATA1':'r'}
     for line in ['28121','28105','28310','28360','29309','29249']:
          plt.text(df.loc[df['Guide ID'] == line][timepoint]+0.1, df.loc[df['Guide ID'] == line]['Plasmid vs T21 enrichment']+.1, line, horizontalalignment='left', size='medium', color=colordict[df.loc[df['Guide ID'] == line]['sgRNA target'].to_list()[0]], weight='regular')
     plt.savefig('Figures/Timepoints/CompareEnrichmentsCRISPRi_labeledGATA1'+timepoint+'.pdf', transparent=True)
     plt.close()
data.to_csv('Figures/Timepoints/CompareEnrichmentsCRISPRi_chrX.csv')


# with Strand colored
data = df.loc[(df['Chromosome_x']=='chrX') & (df['Guide start'] > 48641136-1000) & (df['Guide start'] < 48658755+700+1000) & (df['cutting specificity score'] > 0.2)]
fig, axs = plt.subplots(4,1,figsize=(4,6))
axs = axs.flatten()
i=0
for timepoint in ['Plasmid vs T07 enrichment', 'T07 vs T21 enrichment', 'Plasmid vs T21 enrichment', 'Plasmid vs T21 dCas9 enrichment']:
     sns.scatterplot(x = 'Cut site', y = timepoint , hue='strand', data = data, alpha = 0.5, linewidth = 0, ax=axs[i], marker='.')
     i+=1

plt.savefig('Figures/Timepoints/CRISPRipeakwidth with GS CFD filter.pdf', transparent=True)
plt.close()


fig, axs = plt.subplots(1,2,figsize=(5,5))
axs = axs.flatten()
i=0
for timepoint in ['T07 vs T21 enrichment','Plasmid vs T07 enrichment']:
     sns.relplot(y = 'Plasmid vs T21 enrichment', x = timepoint , data = df,  hue = 'sgRNA target', palette=['tab:gray','k','b','m','y','r'], legend = True, ax=axs[i])
     # for line in ['28121','28105','28310','28360','29309','29249']:
     #      plt.text(df.loc[df['Guide ID'] == line]['Plasmid vs T21 enrichment']+0.1, df.loc[df['Guide ID'] == line][timepoint]+.2, line, horizontalalignment='left', size='medium', color='black', weight='semibold')
     i+=1

plt.savefig('Figures/Timepoints/CompareEnrichmentsCRISPRi_labeled.pdf', transparent=True)
plt.close()

# g = sns.residplot(x = 'Plasmid vs T21 enrichment', y = 'T07 vs T21 enrichment' , data = df, color = 'tab:gray', label = 'Other sgRNA')
# g = sns.residplot(x = 'Plasmid vs T21 enrichment', y = 'T07 vs T21 enrichment' , data = eGATA1_df, color = 'r', label = 'eGATA1 enhancer')
# g = sns.residplot(x = 'Plasmid vs T21 enrichment', y = 'T07 vs T21 enrichment' , data = eHDAC6_df, color = 'g', label = 'eHDAC6 enhancer')
# g = sns.residplot(x = 'Plasmid vs T21 enrichment', y = 'T07 vs T21 enrichment' , data = GATA1_TSS_df, color = 'b', label = 'GATA1 TSS')
# plt.ylabel('Residual enrichment (T07 vs T21)')
# plt.legend()
# plt.savefig('Figures/Timepoints/ResidualEnrichmentsCRISPRi.pdf', transparent=True) ## BUGGY
# plt.close()

g = sns.PairGrid(df[['Plasmid vs T21 enrichment', 'T07 vs T21 enrichment', 'sgRNA target']], hue='sgRNA target') 
g.map_upper(sns.regplot) 
g.map_lower(sns.residplot) 
g.map_diag(plt.hist) 
g.add_legend() 
g.set(alpha=0.5)
plt.savefig('Figures/Timepoints/CompareEnrichmentsCRISPRi_PairGrid.pdf', transparent=True)
plt.close()

# g = sns.PairGrid(df[['Plasmid vs T21 enrichment', 'T07 vs T21 enrichment', 'cutting specificity score']]) 

g = sns.jointplot(x = 'cutting specificity score', y = 'T07 vs T21 enrichment' , data = df, alpha = 0.3, linewidth = 0, space = 0)
plt.savefig('Figures/Timepoints/GS-CFDvsT07timepoint.pdf', transparent=True)
plt.close()

g = sns.jointplot(x = 'cutting specificity score', y = 'Plasmid vs T21 enrichment' , data = df, alpha = 0.3, linewidth = 0, space = 0)
plt.savefig('Figures/Timepoints/GS-CFDvsPlasmid timepoint.pdf', transparent=True)
plt.close()

g = sns.lmplot(x = 'cutting specificity score', y = 'Plasmid vs T21 enrichment' , data = df,  hue = 'sgRNA target', fit_reg=False, palette=['tab:gray','k','b','m','y','r'], legend = True)
for line in ['28121','28105','28310','28360','29309','29249']:
     plt.text(df.loc[df['Guide ID'] == line]['cutting specificity score']+0.1, df.loc[df['Guide ID'] == line]['Plasmid vs T21 enrichment']+.2, line, horizontalalignment='left', size='medium', color='black', weight='semibold')


plt.xlabel('GuideScan CFD specificity score')
plt.savefig('Figures/Timepoints/GS-CFDvsPlasmid timepoint_labeled.pdf', transparent=True)
plt.close()
