'''
####### Author: Jin Woo Oh  ############
###### For CRISPR WG paper  ############
###### Last updated: 04/16/2023#########
########################################

'''

import sys
import numpy as np
def main(argv = sys.argv):
    if(len(argv) <3):
        print("{0} {rep1} {rep2} .. etc")
        sys.exit()

    cnt_vects = []
    id_vects = [] # for checking if the reps are all in same order
    for i in range(1, len(argv)):
        with open(argv[i], 'r') as ifile:
            cnt_vect = [eval(line.split()[3]) for line in ifile]
            cnt_vects.append(cnt_vect)

            ifile.seek(0)
            id_vect = [line.split()[0] for line in ifile] 
            id_vects.append(id_vect)


    avg_cnts = []
    for i in range(0, len(cnt_vects[0])):
        avg_cnt = np.mean([cnt_vects[j][i] for j in range(0, len(cnt_vects))])
        avg_cnts.append(avg_cnt)

        ids = [id_vects[j][i] for j in range(0, len(id_vects))]
        if(len(set(ids)) > 1):
            print("ERROR: input reps in not uniform order")
            sys.exit()
            print(ids)

    with open(argv[1], 'r') as ifile:
        for i, line in enumerate(ifile):
            words = line.split()
            words[3] = str(avg_cnts[i])
            print("\t".join(words))
    

main()
