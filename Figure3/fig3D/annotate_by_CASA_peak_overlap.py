'''
####### Author: Jin Woo Oh: ############
###### For CRISPR WG paper  ############
###### Last updated: 04/16/2023#########
########################################

head GATA_HCRFF_x200.log2FC_0.bedgraph

chr12   54300767        54300767        -0.795300896342266
chr12   54300811        54300811        -0.5280892454653947
chr19   12887237        12887237        -0.5557240249617239
chr3    128487572       128487572       -1.0690192470278501
chr3    128487838       128487838       -0.34549656602576667

head ENCFF413WYU.bed
chrX    48782597        48783497        chrX:48782597-48783497:.        2.231263465     .       chrX:48782597-48783497:.        chrX    4878658948786590 +       NA      NA      GATA1   ENSG00000102145 NA      NA      TRUE    NA      NA      NA      NA      NA      TRUE    NA
chrX    48786097        48786997        chrX:48786097-48786997:.        3.708137016     .       chrX:48786097-48786997:.        chrX    4878658948786590 +       NA      NA      GATA1   ENSG00000102145 NA      NA      TRUE    NA      NA      NA      NA      NA      TRUE    NA
chrX    48800197        48801297        chrX:48800197-48801297:.        2.398983455     .       chrX:48800197-48801297:.        chrX    487865894
'''

import sys

peaks = []
with open(sys.argv[2], 'r') as ifile:
    for line in ifile:
        words = line.split()
        peaks.append([words[0], int(words[1]), int(words[2])])

with open(sys.argv[1]) as ifile:
    for line in ifile:
        words = line.split()
        chr, begin, end = words[0], int(words[1]), int(words[2])

        found = False
        for peak in peaks:
            if(peak[0] == chr and (peak[1] <begin and end < peak[2])):
                print(line.rstrip() + "\t" + "T")
                found = True
                break
        if(not found):
            print(line.rstrip() + "\t" + "F")


