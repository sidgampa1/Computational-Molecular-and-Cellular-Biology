import sys
from Bio import AlignIO
from collections import Counter
import numpy as np
import math
import itertools

## check if correct command line args given
if len(sys.argv) != 2:
    print("Incorrect number of arguments, try again")
    sys.exit()

try:
    alignment = AlignIO.read(sys.argv[1], "stockholm")
    align_array = np.array([list(rec) for rec in alignment], np.character)
except IOError:
    sys.stderr.write("Incorrect type of file, please try again")
    sys.exit()

column_dist = {} ## key = column index, value = dictionary of nuc:probability
column_entropy = {} ## key = column index, value = column entropy
colpair_dist = {} ## key = column pair (i, j), value = dic of nuc pair: entropy
colpair_mi = {} ## key = column pair (i, j), value = mutual information
nucs = ['A', 'U', 'C', 'G', '-']
nuc_pairs = list(itertools.product(nucs, repeat = 2)) ##list of nuc_pairs
seq_columns = align_array.T

## count nuc prob for each column
column_index = 0
for column in seq_columns:
    column = column.tolist()
    prob_dic = Counter(column)
    total = sum(prob_dic.itervalues())
    prob_dic = {k: v/float(total) for k, v in prob_dic.iteritems()}
    for nuc in nucs:
        if nuc not in prob_dic:
            prob_dic[nuc] = 0.0 ## ensures all nucs accounted for
    column_dist[column_index] = prob_dic
    column_index += 1


## calc entropy for each column
for index, column in column_dist.iteritems():
    entropy = 0
    for nuc, prob in column.iteritems():
        if prob > 0:
            entropy += prob * math.log(prob, 2)
    column_entropy[index] = -entropy


## calc column pair probability
nuc_pairs = ["".join(tup) for tup in nuc_pairs]
column_pairs = list(itertools.combinations(range(0, len(column_entropy)), r=2))
colpair_dist = {cols: {} for cols in column_pairs}
nuc_dist_col1 = {nuc: [] for nuc in nucs}
nuc_dist_col2 = {nuc: [] for nuc in nucs}
for cols in column_pairs:
    col_length = len(seq_columns[cols[0]])

    for nuc in nuc_dist_col1:
        nuc_dist_col1[nuc] = [i for i, x in enumerate(seq_columns[cols[0]].tolist()) if x == nuc]
        nuc_dist_col2[nuc] = [i for i, x in enumerate(seq_columns[cols[1]].tolist()) if x == nuc]

    for nucs in nuc_pairs:
        nuc1_dist = nuc_dist_col1[nucs[0]]
        nuc2_dist = nuc_dist_col2[nucs[1]]
        for index in nuc1_dist:
            if index in nuc2_dist:
                if nucs in colpair_dist[cols]:
                    colpair_dist[cols][nucs] += 1
                if nucs not in colpair_dist[cols]:
                    colpair_dist[cols][nucs] = 1

        if nucs not in colpair_dist[cols]:
            colpair_dist[cols][nucs] = 0

        colpair_dist[cols][nucs] /= float(col_length)


## calculate mutual info
for cols, joint_probs in colpair_dist.iteritems():
    nuc_probs_col0 = column_dist[cols[0]]
    nuc_probs_col1 = column_dist[cols[1]]
    col_mi = 0

    for nucs, nucs_prob in joint_probs.iteritems():
        nuc1_prob = nuc_probs_col0[nucs[0]]
        nuc2_prob = nuc_probs_col1[nucs[1]]
        if nuc1_prob > 0 and nuc2_prob > 0 and nucs_prob > 0:
            col_mi += nucs_prob * math.log((nucs_prob / float(nuc1_prob * nuc2_prob)), 2)

    colpair_mi[cols] = col_mi


## round the mutual info and entropy information
for key in colpair_mi.keys():
    colpair_mi[key] = round(colpair_mi[key], 6)

for key in column_entropy.keys():
    column_entropy[key] = round(column_entropy[key], 6)

## sort mutual information and entropy values
colpair_mi_vals = sorted(colpair_mi.values())
colpair_mi_vals.reverse()

column_entropy_vals = sorted(column_entropy.values())

## print out the 10 lowest entropy columns and highest 50 mutual info
for entropy in column_entropy_vals[0:10]:
    index = column_entropy.values().index(entropy)
    column = column_entropy.keys()[index]
    print column
    column_entropy[column] = 0

for mi in colpair_mi_vals[0:50]:
    index = colpair_mi.values().index(mi)
    colpair = colpair_mi.keys()[index]
    print(str(colpair[0]) + "," + str(colpair[1]))
    colpair_mi[colpair] = 0
