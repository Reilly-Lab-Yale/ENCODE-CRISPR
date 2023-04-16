'''
####### Author: Jin Woo Oh  ############
###### For CRISPR WG paper  ############
###### Last updated: 04/16/2023#########
########################################


input file format:
chr10   107843639       107843642       NA|chr10:107843620-107843639:+  355     +       chr10:107843620-107843639:+     chr8    127736231       127736232       +       MYC     ENSG00000136997 GGAATGATCACAGTCTAAAC    GGAATGATCACAGTCTAAAC    negative_control        ST

name    PAM_ID  raw_log2FC      ztransformed_by_neg_control     ztransformed_by_all_guides      normalized_by_mean_count1       normalized_by_mean_count2       normalized_by_tot_count
NA|nt_11        NA:NA-NA:NA     1.4818690077570527      2.2807148116314515      0.7180969205514264      0.14504802392647065     1.1016449190720878     1.1016449190720876

'''
import sys
import numpy as np
def extract_PAM_ID(line):
    words = line.split()
    pam_id = words[0] + ":" + words[1] + "-" + words[2] + ":" + words[5]
    return pam_id
def main(argv = sys.argv):
    if(len(argv) != 4):
        print("Usage: {0} {gRNA quant 1 (e.g. T0)} {gRNa quant 2 (e.g. T14)} {ofile prefix}")
        sys.exit()
    ifile1 = open(argv[1], 'r')
    ifile2 = open(argv[2], 'r')

    prefix = argv[3]

    # Generate N x 3 matrix, N = #gRNA. : [[perturbation ID, type, count]]
    # [['MYC|chr12:54300748-54300767:+', chr12:54300767-54300770:+, 'targeting' 517] ..]
    pID_to_cnt_1 = np.array([[line.split()[3], extract_PAM_ID(line), line.split()[-2], int(line.split()[4])] for line in ifile1], dtype=object)
    pID_to_cnt_2 = np.array([[line.split()[3], extract_PAM_ID(line),  line.split()[-2], int(line.split()[4])] for line in ifile2], dtype=object)

    #if((pID_to_cnt_1[:,0] != pID_to_cnt_2[:,0]).all()):
    if(not np.array_equal(pID_to_cnt_1[:,0], pID_to_cnt_2[:,0])):
        print("gRNA order in file1 and and file2 don't match")
        sys.exit()


    cnt1, cnt2 = pID_to_cnt_1[:,3].astype('float64'), pID_to_cnt_2[:,3].astype('float64')

    # compute log2 (pre-normalization)
    log2FC = np.log2(((cnt1 + 1.0) / (cnt2 + 1.0)))

    # normalize using all
    mu, sig = log2FC.mean(), log2FC.std() 
    log2FC_all_norm = (log2FC - mu)/sig

    # normalize by mean count
    log2FC_mean_normed1 = np.log2(((cnt1/(np.mean(cnt1)) + 1.0) / \
                                  (cnt2/(np.mean(cnt2)) + 1.0)))

    # slightly differnt version suggested by mike. normalize by mean count 
    log2FC_mean_normed2 = np.log2(((cnt1+1.0)/(np.mean(cnt1))) / \
                                  ((cnt2+1.0)/(np.mean(cnt2))))

    # tot normed (suggested by josh)
    log2FC_tot_normed = np.log2(((cnt1+1.0)/(np.sum(cnt1))) / \
                                  ((cnt2+1.0)/(np.sum(cnt2))))

    print("\t".join([str(np.mean(cnt1)), str(np.mean(cnt2))]))
    with open(prefix + ".log2FC_summary", 'w') as ofile:
        ofile.write("\t".join(["name", "PAM_ID","raw_log2FC",  "ztransformed_by_all_guides", "normalized_by_mean_count1", "normalized_by_mean_count2", "normalized_by_tot_count"]) + '\n')
        for i in range(0, len(log2FC)):
            if("NA|.:.-." in pID_to_cnt_1[i,0]):
                    continue
            else:
                ofile.write("\t".join([pID_to_cnt_1[i,0], pID_to_cnt_1[i,1] \
                               ,str(log2FC[i]),  str(log2FC_all_norm[i]), str(log2FC_mean_normed1[i]), str(log2FC_mean_normed2[i]),  str(log2FC_tot_normed[i])]) + '\n')
            
            
    ifile1.close()
    ifile2.close()
    
    

main()
