#!/usr/bin/env python
'''
The script is used to do PAM (Partitioning Around Medoids) clustering.
The script has not been finalized.
The input and output path and the pam clustering number should be replaced manually.
'''

from sys import argv
import pandas as pd
import kmedoids
try:
    import cPickle as pickle
except ModuleNotFoundError:
    import pickle


# load distance matrix
with open("/data/yanw3/GutFunPROJ/decarboxylase/emapper/GroupII/neighbor/data_dist_PY_index_reduce5.pkl", 'rb') as inp:
    df = pickle.load(inp)

# set a clustering number
Knum = 500

c = kmedoids.fasterpam(df, Knum)
print("Loss is:", c.loss)

print(c, file=open(f"/data/yanw3/GutFunPROJ/decarboxylase/emapper/GroupII/neighbor/pam/pam_cluster{Knum}.txt", "w"))


names = pd.DataFrame(df.index)
dfc = pd.DataFrame(c.labels,columns=['cluster'])
dfc['ID'] = df.index
col = ['ID','cluster']
dfc = dfc[col]


dfc.to_csv(f"/data/yanw3/GutFunPROJ/decarboxylase/emapper/GroupII/neighbor/pam/pam_cluster{Knum}_list.tsv", sep='\t',header=False, index=False)
