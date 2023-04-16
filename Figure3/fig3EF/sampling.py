'''
####### Author: Jin Woo Oh  ############
###### For CRISPR WG paper  ############
###### Last updated: 04/16/2023#########
########################################


chrX    48686620        48686623        NA|chrX:48686620-48686623       20      -       chrX:48686620-48686623  NA      NA      NA      NA   NA       NA      GGGTCTATGGAGAAGTGGG     GGGGTCTATGGAGAAGTGGG    targeting       NA
chrX    48686624        48686627        NA|chrX:48686624-48686627       223     -       chrX:48686624-48686627  NA      NA      NA      NA   NA       NA      AGTAGGGTCTATGGAGAAG     GAGTAGGGTCTATGGAGAAG    targeting       NA
chrX    48686632        48686635        NA|chrX:48686632-48686635       184     -       chrX:48686632-48686635  NA      NA      NA      NA   NA       NA      AGGGTTCAAGTAGGGTCTA     GAGGGTTCAAGTAGGGTCTA    targeting       NA
chrX    48686640        48686643        NA|chrX:48686640-48686643       217     -       chrX:48686640-48686643  NA      NA      NA      NA   N


# simulate biorepr correlation at each seq depth
# return average dropout rate, averaged between two bioreps 
'''

import sys
import numpy as np
from scipy import stats
from scipy.stats.stats import pearsonr   


# return n_simul x gRNA matrix, M[i,j] = numberof jth gRNA in i'th simulation
def gRNA_sampling(cnt, N_s, n_simul):
    matrix = []
    pool = []
    n_gRNAs = len(cnt)
    for i in range(0, len(cnt)):
        for j in range(0, int(cnt[i])):
            pool.append(i)

    for i in range(0, n_simul):
   #     np.random.seed(i)
        sampled = np.random.choice(pool, N_s, replace = True)

        # generate count histograms of simulation result
        unique, counts = np.unique(sampled, return_counts=True)
        cdict = dict(zip(unique, counts))


        # convert count dictionary (gRNA index --> count) to an array
        cnts = []
        for i in range(0, n_gRNAs):
            try:
                cnts.append(cdict[i])
            except KeyError: # if not in the dictionary, means count was zero
                cnts.append(0)
        matrix.append(cnts)
    return np.array(matrix)

# given n_simul x gRNA matrix, return a vector of size n_simul to return fraction of guides that are dropped out 
def dropout(matrix1, matrix2, threshold):
    n_simul = len(matrix1)
    n_guide = len(matrix1[0])

    dropout_rates = []
    for i in range(0, n_simul):
         dropout_rates.append(np.sum(  np.logical_or(matrix1[i]<threshold, matrix2[i]<threshold) ) / n_guide)

    return np.array(dropout_rates)


# given simulation matrix ( n_simul x gRNA matrix, M[i,j] = numberof jth gRNA in i'th simulation)
# compute log2FC for each n_simul -> return a single matrix in same dimension
def log2FC(Ma, Mb):
    log2FC_mat = []
    for i in range(0, len(Ma)):
        log2FC = np.log2((1.0 + Ma[i]) / (1.0 + Mb[i]))
        log2FC_mat.append(log2FC)
    return np.array(log2FC_mat)


def main(argv = sys.argv):
    if(len(argv) != 8):
        print("{0} {encode formatted guide quant 1: brep1} {encode formatted guide quant 2: brep1} {encode formatted guide quant 1: brep2} {encode formatted guide quant 2: brep2}  {number of simulations} {dropout threshold}  {ofile prefix}")
        sys.exit()
    sa_cnt_brep1, sb_cnt_brep1 = [], []
    with open(argv[1], 'r') as ifile:
        for line in ifile:
            sa_cnt_brep1.append(int(line.split()[4]))
    with open(argv[2], 'r') as ifile:
        for line in ifile:
            sb_cnt_brep1.append(int(line.split()[4]))

    sa_cnt_brep2, sb_cnt_brep2 = [], []
    with open(argv[3], 'r') as ifile:
        for line in ifile:
            sa_cnt_brep2.append(int(line.split()[4]))
    with open(argv[4], 'r') as ifile:
        for line in ifile:
            sb_cnt_brep2.append(int(line.split()[4]))

    n_sim = int(argv[5])
    dropout_threshold = int(argv[6])
    pref = argv[7]

    sa_cnt_brep1 = np.array(sa_cnt_brep1)
    sb_cnt_brep1 = np.array(sb_cnt_brep1)
    sa_cnt_brep2 = np.array(sa_cnt_brep2)
    sb_cnt_brep2 = np.array(sb_cnt_brep2)


    dropout_filter = np.logical_and(np.logical_and(sa_cnt_brep1>0, sb_cnt_brep1>0), np.logical_and(sa_cnt_brep2>0, sb_cnt_brep2>0))
    sa_cnt_brep1 = sa_cnt_brep1[dropout_filter]
    sb_cnt_brep1 = sb_cnt_brep1[dropout_filter]
    sa_cnt_brep2 = sa_cnt_brep2[dropout_filter]
    sb_cnt_brep2 = sb_cnt_brep2[dropout_filter]


    tot_gRNA = len(sa_cnt_brep1)
#    depth_a = np.sum(sa_cnt_brep1) / float(tot_gRNA)
#    depth_b = np.sum(sb_cnt_1) / float(tot_gRNA)
#    print("Depth of Samp A: " + str(depth_a))
#    print("Depth of Samp B: " + str(depth_b))
    # turn into sample discrete probability distribution


    #############################  sampling ###########################################
    ofile_corrs = open(pref + ".samp_corrs", 'w', buffering=1)
    ofile_corr_stat = open(pref + ".samp_corrs_stat", 'w', buffering=1)
    ofile_dropout_stat = open(pref + ".samp_dropout_stat", 'w', buffering=1)

    for depth in list(range(1, 1020, 30)) + [5000]:
        N_sample = depth * tot_gRNA
        sim_sa_cnt_brep1 = gRNA_sampling(sa_cnt_brep1,  int(N_sample), n_sim)
        sim_sb_cnt_brep1 = gRNA_sampling(sb_cnt_brep1,  int(N_sample), n_sim)
        sim_sa_cnt_brep2 = gRNA_sampling(sa_cnt_brep2,  int(N_sample), n_sim)
        sim_sb_cnt_brep2 = gRNA_sampling(sb_cnt_brep2,  int(N_sample), n_sim)

        d1 = dropout(sim_sa_cnt_brep1, sim_sb_cnt_brep1, dropout_threshold) # vector of length n_sim
        d2 = dropout(sim_sa_cnt_brep2, sim_sb_cnt_brep2, dropout_threshold)
        dropout_rate = (d1+d2)/2

        log2FCs_biorep1 = log2FC(sim_sa_cnt_brep1, sim_sb_cnt_brep1)
        log2FCs_biorep2 = log2FC(sim_sa_cnt_brep2, sim_sb_cnt_brep2)

        corrs = []
        for i in range(0, len(log2FCs_biorep1)):
            for j in range(0, len(log2FCs_biorep1)):
                corrs.append(pearsonr(log2FCs_biorep1[i], log2FCs_biorep2[j])[0])
        
             
        ofile_corrs.write(str(depth) + "\t".join([str(corr) for corr in corrs]) + '\n')
        ofile_corr_stat.write("\t".join([str(depth), str(np.mean(corrs)), str(np.std(corrs)), str(min(corrs)), str(max(corrs))]) + '\n')
        ofile_dropout_stat.write("\t".join([str(depth), str(np.mean(dropout_rate)), str(np.std(dropout_rate)), str(min(dropout_rate)), str(max(dropout_rate))]) + '\n')


    ofile_corrs.close()
    ofile_corr_stat.close()
    ofile_dropout_stat.close()
main()

