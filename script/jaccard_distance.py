#!/usr/bin/env python
'''
The script is used to calculate the jaccard distance between two neighborhood region based on the overlap of eggnod IDs.
The script has not been finalized. The input and out path should be replaced manually.
'''

from itertools import product, starmap
import pandas as pd

def jaccard_similarity(list1, list2):
    intersection = len(set(list1).intersection(list2))
    union = (len(set(list1)) + len(set(list2))) - intersection
    return 1- intersection / union

my_combine_neibr0_reduce4 = pd.read_csv("/data/yanw3/GutFunPROJ/decarboxylase/emapper/GroupII/neighbor/combine.neibr0_reduce5.txt",sep='\t',index_col=False,header=None)
my_list = pd.read_csv("/data/yanw3/GutFunPROJ/decarboxylase/emapper/GroupII/neighbor/nr_cor_list3.keep.list",sep='\t',index_col=False,header=None)


names = list(my_list.iloc[:,0])
myl = []
for i in list(range(0,len(names))):
    tmp = list(my_combine_neibr0_reduce4.loc[my_combine_neibr0_reduce4[3] == names[i],6])
    myl.append(tmp) 


inputs2 = product(myl, myl)
result2 = list(starmap(jaccard_similarity, inputs2))

def group_full(it):
    iterators = [iter(it)] * 31527
    return zip(*iterators)


ca = list(group_full(result2))
daDist = pd.DataFrame(ca)
daDist.columns = names
daDist.index = names

# inspect
daDist.iloc[1:5, 1:5]

# save object
try:
    import cPickle as pickle
except ModuleNotFoundError:
    import pickle

def save_object(obj, filename):
    with open(filename, 'wb') as outp:  # Overwrites any existing file.
        pickle.dump(obj, outp, pickle.HIGHEST_PROTOCOL)

save_object(daDist, '/data/yanw3/GutFunPROJ/decarboxylase/emapper/GroupII/neighbor/data_dist_PY_index_reduce5.pkl')

# write out the dist data
daDist.to_csv(f"/data/yanw3/GutFunPROJ/decarboxylase/emapper/GroupII/neighbor/data_dist_PY_index_reduce5.tsv", header=True, index=True,sep="\t")
