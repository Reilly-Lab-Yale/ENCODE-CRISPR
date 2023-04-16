
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



# file index for best isoforms (pol2 signal). For duplicate lines, choose any one (alternative spliced with same TSS TES).

(base) [joh27@dolomite using_pol2_chip_TSS_refgene]$ head K562_best_isoform_indices.txt
COG6    7257
NT5DC4  3477
SREBF2-AS1      86374
MIR5003 54830
FGF12   6280
FGF12   33557
TMEM161B-AS1    53273
TMEM161B-AS1    53274





(base) [joh27@troctolite strand_analysis]$ head ENCFF274YGF.bed
chr1    181400  181530  .       0       .       0.299874        -1      -1      75
chr1    778660  778800  .       0       .       14.1383 -1      -1      75
chr1    779137  779200  .       0       .       0.33144 -1      -1      75
chr1    827460  827554  .       0       .       3.38384 -1      -1      75
chr1    842880  843060  .       0       .       1.64457 -1      -1      75


^ narrow peaks. extend from center by 1k bp

'''
import sys
def main(argv = sys.argv):
    if(len(argv) != 10):
        print("{0} {bedgraph} {.gtf} {gene name} {K562_best_isoform_indices.txt} {K562 DHS bed file} {ofname 1: gene body} {ofname 2: promoter} {ofile 3: null regions} {ofile 4: TSS flank}")
        sys.exit()

    bfile = open(argv[1], 'r')
    gfile = open(argv[2], 'r')
    gene = argv[3]
    gene_to_index = {line.split('\t')[0] : int(line.split('\t')[1]) for line in open(argv[4])}
    dhs_bfile = open(argv[5], 'r')
    ofile1 = open(argv[6], 'w')
    ofile2 = open(argv[7], 'w')
    ofile3 = open(argv[8], 'w')
    ofile4 = open(argv[9], 'w')
    # K562 DHSs
    chrom_to_dhs = dict()
    for line in dhs_bfile:
        words = line.split()
        chrom, begin, end = words[0], eval(words[1]), eval(words[2])
        center = (begin+end)//2
        begin, end = center - 500, center + 500

        try:
            chrom_to_dhs[chrom].append([begin,end])
        except KeyError:
            chrom_to_dhs[chrom] = [[begin,end]]

    gene_range = []
    promoter_range = []
    gfile_lines = gfile.readlines()
    words = gfile_lines[gene_to_index[gene]].split('\t')
    gchrom, start, end = words[2], int(words[4]), int(words[5])
    true_gene_range = [start, end]

    strand = words[3]
    if(strand == '+'):
        TSS = start
        gene_range = [start + 2000, end] # more stringent gene range to separate out promoter signal
        promoter_range = [start - 2000, start]
    elif(strand == '-'):
        TSS = end
        gene_range = [start, end - 2000]
        promoter_range = [end, end + 2000]
    else:
        print(strand)
        print("??")
        sys.exit()
    print(gene_range)
    print(promoter_range)
    for line in bfile:
        words = line.split()
        chrom, coord = words[0], int(words[1])
        in_gene_body = gchrom == chrom and coord > gene_range[0] and coord < gene_range[1]
        in_true_gene_body = gchrom == chrom and coord > true_gene_range[0] and coord < true_gene_range[1] 
        in_promoter = gchrom == chrom and coord > promoter_range[0] and coord < promoter_range[1]

        in_flank = gchrom == chrom and coord > TSS-5000 and coord < TSS + 5000

        if(in_gene_body):
            ofile1.write(line)

        elif(in_promoter):
            ofile2.write(line)
        
        else: # outside gene body and promoter 
            if(not in_true_gene_body and chrom == gchrom and abs(TSS - coord) < 500000):
                in_dhs = False
                dhs_list = chrom_to_dhs[chrom]
                for dhs in dhs_list:
                    if(coord>dhs[0] and coord<dhs[1]):
                        in_dhs = True
                        break
                if(not in_dhs):
                    ofile3.write(line)

        if(in_flank):
            ofile4.write(line)

    ofile1.close()
    ofile2.close()
    ofile3.close()
    ofile4.close()
    dhs_bfile.close()
    bfile.close()
    gfile.close()

main()
