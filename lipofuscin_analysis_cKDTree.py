import os, sys
import numpy as np
import pandas as pd
#from glob import glob 
#from skimage.measure import regionprops
from skimage.io import imread, imsave
from os.path import abspath, dirname
from scipy.spatial import distance, cKDTree

spot_c0_path      = sys.argv[1]
spot_c1_path      = sys.argv[2]
outdir            = sys.argv[3]
round             = sys.argv[4]          

spot_c0 = np.loadtxt(spot_c0_path)
spot_c1 = np.loadtxt(spot_c1_path,delimiter=',')

vox=[0.23,0.23,0.42]
c0=spot_c0[:,:3]/vox
c1=spot_c1[:,:3]/vox

neighbor_radius   = 3
kdtree_c0 = cKDTree(c0)
kdtree_c1 = cKDTree(c1)
neighbors = kdtree_c0.query_ball_tree(kdtree_c1, neighbor_radius)

no_neighbors = 0
one_neighbor = 0
more_neighbors = 0
max_neighbors = 0
for nnn in neighbors:
    if len(nnn) == 0: no_neighbors += 1
    if len(nnn) == 1: one_neighbor += 1
    if len(nnn) > 1 : more_neighbors += 1
    if len(nnn) > max_neighbors: max_neighbors = len(nnn)

print(no_neighbors, one_neighbor, more_neighbors, max_neighbors)

neighbors_num = np.array([len(x) for x in neighbors]).sum()
pA = np.empty((neighbors_num, 3))
pB = np.empty((neighbors_num, 3))
pAind = np.empty(neighbors_num, dtype=np.uint32)
pBind = np.empty(neighbors_num, dtype=np.uint32)

p_ind = 0
for c0_ind, nnn in enumerate(neighbors):
    if len(nnn) == 0:
        continue
    for c1_ind in nnn:
        pA[p_ind]     = c0[c0_ind, :3]
        pB[p_ind]     = c1[c1_ind, :3]
        pAind[p_ind]  = c0_ind
        pBind[p_ind]  = c1_ind
        p_ind += 1

lipo_c0 = spot_c0[pAind]
lipo_c1 = spot_c1[pBind]

true_pos_c0 = np.delete(spot_c0, pAind, axis=0)
true_pos_c1 = np.delete(spot_c1, pBind, axis=0)

np.savetxt(outdir+'/{}_c0_lipo.txt'.format(round), lipo_c0)
np.savetxt(outdir+'/{}_c1_lipo.txt'.format(round), lipo_c1)
np.savetxt(outdir+'/{}_c0_true.txt'.format(round), true_pos_c0)
np.savetxt(outdir+'/{}_c1_true.txt'.format(round), true_pos_c1)
