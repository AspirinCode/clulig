#!/usr/bin/env python

#############
#
# Cluster molecules in MOL files by similarity using Tanimoto coefficients
# based on Open Babel FP2 fingerprints using the DBSCAN algorithm.
#
# 1. Read all MOL files in the current directory
# 2. Calculate all FP2 fingerprints with Open Babel's pybel module
# 3. Construct distance matrix with each element as 1 - similarity
# 4. Cluster using scikit-learn's DBSCAN module
# 5. Print clusters 
#
#############

import sys
from glob import iglob
from sklearn.cluster import DBSCAN
import numpy as np
import pybel as pb

#We want members of a cluster to have over sim_cutoff similarity.
sim_cutoff=0.85
#Number of neighbors required for definition of a core point in a cluster
min_samples=2

if len(sys.argv) == 2:
    try:
        temp = float(sys.argv[1])
        if temp > 0 and temp < 1.0:
            sim_cutoff = temp
    except Exception:
        pass

mols = []
for molfilename in iglob('*.mol'):
    try:
        m = pb.readfile('mol',molfilename).next()
    except Exception:
        print 'Error reading '+molfilename+' Will skip.'
        continue
    try:
        fp = m.calcfp()
    except Exception:
        print 'Error calculating fingerprint for '+molfilename+'. Will skip.'
        continue
    mols.append((m,fp,molfilename))

nmol = len(mols)
if nmol == 0:
    print 'No MOL files found!'
    sys.exit()
else:    
    print str(nmol)+' MOL files found.'

print 'Will cluster with cutoff of '+str(sim_cutoff)

#Distance is defined as 1-Tanimoto coefficient.
#This is the simplest option; it could be optimized.
distance_matrix = np.empty([nmol,nmol])
for i in xrange(nmol):
    for j in xrange(i,nmol):
        distance_matrix[i][j] = 1.0 - (mols[i][1] | mols[j][1])
        distance_matrix[j][i] = distance_matrix[i][j]

clusterer = DBSCAN(eps=(1.-sim_cutoff),metric='precomputed')
clusterer.fit(distance_matrix)

#labels is a list of length nmol with each element as a cluster label.
#outliers not belonging to a cluster are labeled -1
labels = clusterer.labels_
numclusters = len(set(labels)) - (1 if -1 in labels else 0)

print str(numclusters)+' clusters found.'

for i in xrange(numclusters):
    print 'CLUSTER '+str(i)+':'+' ('+str(sum(labels==i))+' members)'
    for j,k in enumerate(labels):
        if k==i:
            print mols[j][2]+' '+mols[j][0].title+' cluster'+str(i)
    print ''

