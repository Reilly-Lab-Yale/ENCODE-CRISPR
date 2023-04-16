'''
A_Bassik
CAPRIN1_HCRFF
CAT_HCRFF
CD164_HCRFF
D_Bassik
...

ERP29_HCRFF_coding_gene_body_ravg.log2FC_3.bedgraph
ERP29_HCRFF_template_gene_body_ravg.log2FC_3.bedgraph
ERP29_HCRFF_coding_promoter_ravg.log2FC_3.bedgraph
ERP29_HCRFF_template_promoter_ravg.log2FC_3.bedgraph


output
samp    coding_gene_body    template_gene_body  coding_prommoter    template_promoter
A_Bassik    0.1  0.4  0.2  0.2 
'''


import sys
import numpy as np
def fname_to_mean(fname):
    vals = []
    with open(fname, 'r') as ifile:
        for line in ifile:
            val = eval(line.split()[-1])
            vals.append(val)
    if(len(vals) > 0):
        return len(vals), np.mean(vals)
    else:
        return len(vals), 0

def main(argv = sys.argv):
    if(len(argv) != 2):
        print("{0} {samp list}")
        sys.exit()

    ifile = open(argv[1], 'r')
    print("\t".join(["samp", "perturbation_type", "readout", "coding_gene_body", "N_cgb", "template_gene_body", "N_tgb", "coding_promoter", "N_cp", "template_promoter", "N_tp"]))
    for line in ifile:
        samp,type, readout = line.rstrip().split("\t")
        ncg, cg = fname_to_mean(samp + "_coding_gene_body_ravg.log2FC_Z.bedgraph")
        ntg, tg = fname_to_mean(samp + "_template_gene_body_ravg.log2FC_Z.bedgraph")
        ncp, cp = fname_to_mean(samp + "_coding_promoter_ravg.log2FC_Z.bedgraph")
        ntp, tp = fname_to_mean(samp + "_template_promoter_ravg.log2FC_Z.bedgraph")
        print("\t".join([samp, type, readout, str(cg), str(ncg), str(tg), str(ntg), str(cp), str(ncp), str(tp), str(ntp)]))
    ifile.close()
main()
