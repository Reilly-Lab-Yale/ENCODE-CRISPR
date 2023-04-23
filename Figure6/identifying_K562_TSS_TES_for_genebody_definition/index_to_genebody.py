'''
GENE_A  index


[joh27@dolomite K562_TSS_for_kuei]$ head /mnt/t/data0/joh27/projects/ENCODE/ref_genome/gene_annotation/refGene.txt
585     NR_024540       chr1    -       14361   29370   29370   29370   11      14361,14969,15795,16606,16857,17232,17605,17914,18267,24737,29320,       14829,15038,15947,16765,17055,17368,17742,18061,18366,24891,29370,      0       WASH7P  unk     unk     -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
585     NR_106918       chr1    -       17368   17436   17436   17436   1       17368,  17436,  0       MIR6859-1       unk     unk     -1,
585     NR_107062       chr1    -       17368   17436   17436   17436   1       17368,  17436,  0       MIR6859-2       unk     unk     -1,

gfile_lines = gfile.readlines()
words = gfile_lines[gene_to_index[gene]].split('\t')
gchrom, start, end = words[2], int(words[4]), int(words[5])
'''

import sys
argv = sys.argv

gtf = open(argv[2], 'r').readlines()
genes = []
with open(argv[2], 'r') as ifile:
    for line in ifile:
        words = line.split()
        genes.append([words[12], words[2], words[3], words[4], words[5]])
print("\t".join(["gene_name","chr","strand","gene_start","gene_end"]))
with open(argv[1], 'r') as ifile:
    for line in ifile:
        words = line.rstrip().split('\t')
        gene, index = words[0], int(words[1])
        ginfo = genes[index]
        if(ginfo[0] != gene):
            print("mismatch")
            print(ginfo)
            print(gene)
            sys.exit()

        print("\t".join(ginfo))
