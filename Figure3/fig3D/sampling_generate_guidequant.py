'''
####### Author: Jin Woo Oh  ############
###### For CRISPR WG paper  ############
###### Last updated: 04/16/2023#########
########################################


ifile format 
chrX    48686620        48686623        NA|chrX:48686620-48686623       20      -       chrX:48686620-48686623  NA      NA      NA      NA   NA       NA      GGGTCTATGGAGAAGTGGG     GGGGTCTATGGAGAAGTGGG    targeting       NA
chrX    48686624        48686627        NA|chrX:48686624-48686627       223     -       chrX:48686624-48686627  NA      NA      NA      NA   NA       NA      AGTAGGGTCTATGGAGAAG     GAGTAGGGTCTATGGAGAAG    targeting       NA
chrX    48686632        48686635        NA|chrX:48686632-48686635       184     -       chrX:48686632-48686635  NA      NA      NA      NA   NA       NA      AGGGTTCAAGTAGGGTCTA     GAGGGTTCAAGTAGGGTCTA    targeting       NA
chrX    48686640        48686643        NA|chrX:48686640-48686643       217     -       chrX:48686640-48686643  NA      NA      NA      NA   N


'''

import sys
import numpy as np
from scipy import stats
from scipy.stats.stats import pearsonr   


# return vector of gRNA cnt
def gRNA_sampling(cnt, N_s, rseed):
    matrix = []
    pool = []
    n_gRNAs = len(cnt)
    for i in range(0, len(cnt)):
        for j in range(0, int(cnt[i])):
            pool.append(i)

    np.random.seed(rseed)
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
    return np.array(cnts)

def main(argv = sys.argv):
    if(len(argv) != 5):
        print("{0} {encode formatted guide quant}  {target sequencing depth} {number of simulations}   {ofile prefix}")
        sys.exit()

    input_lines = []
    s_cnt = []
    with open(argv[1], 'r') as ifile:
        for line in ifile:
            input_lines.append(line.rstrip().split('\t'))
            s_cnt.append(int(line.split()[4]))

    depth = int(argv[2])
    n_sim = int(argv[3])
    pref = argv[4]

    s_cnt = np.array(s_cnt)




    tot_gRNA = len(s_cnt)

    
    N_sample = depth * tot_gRNA
    print(tot_gRNA)
    for i in range(0, n_sim):
        sim_s_cnt = gRNA_sampling(s_cnt,  int(N_sample), i)

        of = open(pref + "_sim_" + str(i+1) + "_guidequant.txt", 'w')

                
        for j in range(0, tot_gRNA):
            output_line = input_lines[j][:]
            output_line[4] = str(sim_s_cnt[j])

            of.write("\t".join(output_line) + "\n")
        of.close()



main()

