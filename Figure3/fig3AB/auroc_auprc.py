'''
####### Author: Jin Woo Oh  ############
###### For CRISPR WG paper  ############
###### Last updated: 04/16/2023#########
########################################

# print out AUROC, AUPRC
# save PRC curve

given a file like this :

ifile format 

chr12   54300767        54300767        0.363057970847217       F
chr12   54300811        54300811        0.015926877937254433    F
chr12   54300999        54300999        0.49195603201099186     F
chr12   54301042        54301042        -1.117941257371417      F


compute auroc and auprc.

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

    if(len(argv) != 4):
        print("{0} {ifile} {num intervals} {ofile}")
        sys.exit()

    vals = []; labels = []
    N = int(argv[2])

    with open(argv[1], 'r') as ifile:
        for line in ifile:
            words = line.split()
            vals.append(eval(words[3]))
            labels.append(words[4])

    vals, labels = np.array(vals), np.array(labels)

    inf, sup = np.min(vals) , np.max(vals) - 0.0001 

    precisions = []
    recalls = []
    specificities = []

    for t in np.linspace(sup, inf, N):
        TP = sum(np.logical_and(vals > t, labels == "T"))
        TN = sum(np.logical_and(vals < t, labels == "F"))
        FP = sum(np.logical_and(vals > t, labels == "F"))
        FN = sum(np.logical_and(vals < t, labels == "T"))
        
    
        precisions.append(TP / (TP + FP))
        recalls.append(TP / (TP + FN))
        specificities.append(TN / (TN + FP))

    precisions, recalls, specificities = np.array(precisions), np.array(recalls), np.array(specificities)  
    print(AUC(recalls, precisions))
    print(AUC(1 - specificities, recalls)) 

    m = len(precisions)
    print(len(recalls))
    print(len(specificities))
    print(len(precisions))
    with open(argv[-1], 'w') as ofile:
        for i in range(0,m):
            r = str(recalls[i])
            p = str(precisions[i])
            ofile.write("\t".join([p,r]) + '\n')
main()


