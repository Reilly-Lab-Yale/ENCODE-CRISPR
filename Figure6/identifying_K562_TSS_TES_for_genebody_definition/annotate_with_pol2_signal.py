'''
(base) [joh27@dolomite K562_TSS_for_kuei]$ head /mnt/t/data0/joh27/projects/ENCODE/FCC/crispr_group/final_ENCODE_formatted_data/log2FC/strand_analysis/isoform_index_with_K562_pol2_TSS.out
AADAC   528     0.039281
MIR4660 45158   1.1914
PPWD1   67973   288.286
PPWD1   67977   288.286
KSR2    39693   0.00173516

(base) [joh27@dolomite K562_TSS_for_kuei]$ head K562_TSS_selected_by_higest_Pol2_TSS.out
gene_name       chr     TSS
AADAC   chr3    151814115
MIR4660 chr8    9048444
PPWD1   chr5    65563238
KSR2    chr12   117968990

'''

import sys

gene_to_sig = dict()
with open("/mnt/data0/joh27/projects/ENCODE/FCC/crispr_group/final_ENCODE_formatted_data/log2FC/strand_analysis/isoform_index_with_K562_pol2_TSS.out", 'r') as ifile:
    for line in ifile:
        words = line.rstrip().split('\t')
        gene_to_sig[words[0]] = words[-1]

print("\t".join(["gene_name","chr","TSS","TSS_POL2_sig"]))
with open('K562_TSS_selected_by_higest_Pol2_TSS.out', 'r') as ifile:
    for line in ifile:
        words = line.rstrip().split('\t')
        sig = gene_to_sig[words[0]]
        print("\t".join(words + [sig]))

