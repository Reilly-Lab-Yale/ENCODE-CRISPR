'''

####### Author: Jin Woo Oh  ############
###### For CRISPR WG paper  ############
###### Last updated: 04/16/2023#########
########################################
:wq
Purpose:
Input file has three columns for different types of log2fc. Choose one column and make bedgraph.


input file format:
ame    PAM_ID  raw_log2FC      ztransformed_by_neg_control     ztransformed_by_all_guides      normalized_by_mean_count
GATA1|chr12:54300749-54300767:+ chr12:54300767-54300770:+       -0.6898171127857369     -1.3909202490331116     -1.1558450982161879     -0.5405896621752653
GATA1|chr12:54300793-54300811:+ chr12:54300811-54300814:+       0.14274017211608214     -0.3609625301812802     -0.3592771390999452     -0.15369424933173065
GATA1|chr12:54300981-54300999:+ chr12:54300999-54301002:+       0.022809882908085926    -0.5093284496608929     -0.47402314810507945    -0.21028124214200453



output file format:
chr8    128391200       128391200       1.254248711277702
chr8    126891867       126891867       0.9836344730213816


For output pam coord, use first base pair (if - strand, last coord)
If there are mutliple pam with the same first bp (e.g. coming from pos and neg), average log2FC.
'''
import sys
import numpy as np

def main(argv = sys.argv):
    if(len(argv) != 4):
        print("Usage: {0} {log2FC summary} {log2fc type: e.g. normalized_by_mean_count} {ofile name}")
        sys.exit()

    ifile = open(argv[1], 'r')
    ntype = argv[2]
    header = ifile.readline().split()
    cindex = header.index(ntype)
    ofile = open(argv[3], 'w')

    ploc_to_log2FC = dict()
    for line in ifile:
        words = line.split()
        pID, log2FC = words[1], eval(words[cindex])
        chrom, brange, strand = pID.split(":")

        if(len(chrom.split('_')) > 1): # include only autosome and sex
            continue

        start, end = brange.split("-")

        if(start == "."): # ignore non-coordinates 
            continue 
        if(strand == "+"):
            pcoord = start
        elif(strand == "-"):
            pcoord = str(int(end) - 1)
        else:
            continue
        ploc = chrom + "_" + pcoord
        try:
            ploc_to_log2FC[ploc].append(log2FC)
        except KeyError:
            ploc_to_log2FC[ploc] = [log2FC]

    for ploc in list(ploc_to_log2FC):
        log2FC = np.mean(ploc_to_log2FC[ploc])
        chrom, pcoord = ploc.split("_")
        ofile.write("\t".join([chrom, pcoord, pcoord, str(log2FC)]) + '\n')
    ifile.close()


main()
