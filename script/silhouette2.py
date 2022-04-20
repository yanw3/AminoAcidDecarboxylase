#!/usr/bin/env python
'''
This script is used to calculate the silhouette values to evaluate the clustering.
'''

import pandas as pd
import kmedoids
from sys import argv
import numpy as np

try:
    import cPickle as pickle
except ModuleNotFoundError:
    import pickle

# input the cluster numbers, path to distance matrix, output path, and output prefix.
script, Knum, distFile, outPath, preName = argv
Knum = int(Knum)

# load distance matrix
with open(distFile, 'rb') as inp:
    df = pickle.load(inp)

# df = pd.read_csv(distFile,sep='\t',index_col=0,header=(0))
# calculate k-medoides and silhouette evaluation
c = kmedoids.fasterpam(df, Knum,random_state=19)
silut = kmedoids.silhouette(df, c.labels, samples=False, n_cpu=- 1)

# output
myar = {'n_clusters':[Knum],'silhouette':[silut]}
dac = pd.DataFrame(myar)
dac.to_csv(f'{outPath}/{preName}_{Knum}.tsv', sep='\t',header=True, index=False)
