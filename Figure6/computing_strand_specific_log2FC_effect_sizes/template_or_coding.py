'''
####### Author: Jin Woo Oh  ############
###### For CRISPR WG paper  ############
###### Last updated: 04/16/2023#########
########################################


refGene.txt format
1056    NM_013402       chr11   -       61799628        61817003        61802410        61816929        12      61799628,61802800,61803031,61803362,61803669,61804684,61806663,61810750,61810973,61812470,61813242,61816554,    61802462,61802926,61803111,61803459,61803767,61804761,61806724,61810879,61811075,61812668,61813353,61817003,    0       FADS1   cmpl    cmpl    2,2,0,2,0,1,0,0,0,0,0,0,

(
string  geneName;           "Name of gene as it appears in Genome Browser."
string  name;               "Name of gene"
string  chrom;              "Chromosome name"
char[1] strand;             "+ or - for strand"
uint    txStart;            "Transcription start position"
uint    txEnd;              "Transcription end position"
uint    cdsStart;           "Coding region start"
uint    cdsEnd;             "Coding region end"
uint    exonCount;          "Number of exons"
uint[exonCount] exonStarts; "Exon start positions"
uint[exonCount] exonEnds;   "Exon end positions"
)

'''

import sys
import os 
def main(argv = sys.argv):
    if(len(argv) != 4):
        print("{0} {bedgraph} {gene name} {gene.gtf}")
        sys.exit()

    fname, gene = argv[1], argv[2]
    gfile = open(argv[3], 'r')

    if("pstrand" in fname):
        strand = '+'
    elif("nstrand" in fname):
        strand = '-'
    else:
        print("fname must contain pstrand or nstrand")
        sys.exit()

    gene_strand = '.'
    for line in gfile:
        if(gene in line.split()):
            gene_strand = line.split()[3]
            break
    gfile.close()

    if(gene_strand == '.'):
        print("?")
        sys.exit()

    label = ""
    if(strand == gene_strand):
        label = "template"
    else:
        label = "coding" 

    if("pstrand" in fname):
        newname = fname.replace("pstrand", label)
    else:
        newname = fname.replace("nstrand", label)
    print(fname)
    print(newname)
    print(label)
    os.system("cp " + fname + " " + newname)

main()
