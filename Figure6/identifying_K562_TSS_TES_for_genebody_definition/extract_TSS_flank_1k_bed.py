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
def main(argv = sys.argv):
    if(len(argv) != 2):
        print("{0}  {.gtf}")
        sys.exit()

    gfile = open(argv[1], 'r')

    gene_range = []
    TSSflank_range = []
    for i, line in enumerate(gfile):
        words = line.split()
        start, end = int(words[4]), int(words[5])
        chrom, strand = words[2], words[3]
        gene = words[12]
        if(strand == '+'):
            TSSflank_range = [start - 500, start + 500]
        elif(strand == '-'):
            TSSflank_range = [end-500, end + 500]
        else:
            print(strand)
            print("??")
            sys.exit()

        start = str(TSSflank_range[0])
        end = str(TSSflank_range[1])
        print("\t".join([chrom, start, end , str(i) + ":" + gene]))

main()
