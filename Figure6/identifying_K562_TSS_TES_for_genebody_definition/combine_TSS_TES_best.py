'''
####### Author: Jin Woo Oh  ############
###### For CRISPR WG paper  ############
###### Last updated: 04/16/2023#########
########################################

'''


import sys

def main(argv = sys.argv):
    if(len(argv) != 3):
        print("{0} {isoform_index_with_K562_pol2_TSS.out} {isoform_index_with_K562_pol2_TES.out}")
        sys.exit()


    gene_to_isoforms_highest_TES = dict()
    with open(argv[2], 'r') as ifile:
        for line in ifile:
            words = line.rstrip().split('\t')
            try:
                gene_to_isoforms_highest_TES[words[0]].append(words[1])
            except KeyError:
                gene_to_isoforms_highest_TES[words[0]] = [words[1]]

    gene_to_best_isoforms = dict()
    with open(argv[1], 'r') as ifile:
        for line in ifile:
            words = line.rstrip().split('\t')
            if(words[1] in gene_to_isoforms_highest_TES[words[0]]):
                try:
                    gene_to_best_isoforms[words[0]].append(words[1])
                except KeyError:
                    gene_to_best_isoforms[words[0]] = [words[1]]

    for gene in list(gene_to_best_isoforms):
        for isoform in gene_to_best_isoforms[gene]:
            print("\t".join([gene, isoform]))
main()
