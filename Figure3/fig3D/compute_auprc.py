'''
####### Author: Jin Woo Oh: ############
###### For CRISPR WG paper  ############
###### Last updated: 04/16/2023#########
########################################




# compute AUPRC for table

# output format
100     1000    10      0.8102106596203971

# input format 
chr12   54300767        54300767        0.363057970847217       F
chr12   54300811        54300811        0.015926877937254433    F
chr12   54300999        54300999        0.49195603201099186     F
chr12   54301042        54301042        -1.117941257371417      F


compute  auprc.

'''
import sys
import numpy as np

def AUC(x, y): # trapezoid

    auc = 0
    if(len(x) != len(y)):
        print("vector size mismatch")
        sys.exit()

    dX = [x[i+1] - x[i] for i in range(0, len(x)- 1)]
    for i in range(0, len(y)-1):
        auc += (y[i] + y[i+1]) / 2 * dX[i] 
    return auc
def main(argv = sys.argv):

    if(len(argv) != 7):
        print("{0} {ifile} {num intervals} {cell cov} {seq depth} {sim id} {ofile}")
        sys.exit()

    vals = []; labels = []
    N = int(argv[2])

    cc, sd, sim_id = argv[3], argv[4], argv[5]

    with open(argv[1], 'r') as ifile:
        for line in ifile:
            words = line.split()
            vals.append(eval(words[3]))
            labels.append(words[4])

    vals, labels = np.array(vals), np.array(labels)

    inf, sup = np.min(vals) , np.max(vals) - 0.0001 

    precisions = []
    recalls = []

    for t in np.linspace(sup, inf, N):
        TP = sum(np.logical_and(vals > t, labels == "T"))
        TN = sum(np.logical_and(vals < t, labels == "F"))
        FP = sum(np.logical_and(vals > t, labels == "F"))
        FN = sum(np.logical_and(vals < t, labels == "T"))
        
    
        precisions.append(TP / (TP + FP))
        recalls.append(TP / (TP + FN))

    precisions, recalls = np.array(precisions), np.array(recalls)  
    auprc = AUC(recalls, precisions)
    with open(argv[-1], 'w') as ofile:
        ofile.write("\t".join([cc, sd, sim_id, str(auprc)]))

main()


