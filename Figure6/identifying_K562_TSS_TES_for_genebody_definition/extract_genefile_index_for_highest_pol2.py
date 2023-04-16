'''
####### Author: Jin Woo Oh  ############
###### For CRISPR WG paper  ############
###### Last updated: 04/16/2023#########
########################################

'''

import sys

def main(argv = sys.argv):
    if(len(argv) != 2):
        print("{0} {all_isoform_TSS_upstream_1kb_Pol2_chip.out} ")
        sys.exit()

    id_to_chip = dict()
    with open(argv[1], 'r') as ifile:
        for line in ifile:
            words = line.split('\t')
            id_to_chip[words[0]] = words[-2]
    df = []
    genes = set()
    with open(argv[1] , 'r') as ifile:
        for line in ifile:
            words = line.rstrip().split('\t')
            ID = words[0]
            chip = words[-2]
            gene = ID.split(':')[-1]
            genes.add(gene)
            gene_file_index = ID.split(':')[0]

            df.append([gene, gene_file_index, chip])

    for gene in list(genes):
        subdf = [item for item in df if item[0] == gene]
        max_index = "?"
        max_val = -1
        for i, item in enumerate(subdf):
            if(eval(item[-1]) > max_val):
                max_index = i
                max_val = eval(item[-1])
        if(max_index == "?"):
            print("debug")
            sys.exit()


        # if mutliple TSS with max signal exist, return all
        subdf_max = [item for item in subdf if eval(item[-1]) == max_val]
        for item in subdf_max:
            print("\t".join([gene, item[1], item[2]]))

main()
