#!/usr/bin/env python
# coding: utf-8

# In[2]:


import os, sys
import numpy as np
import pandas as pd
from glob import glob 
from skimage.measure import regionprops
from skimage.io import imread, imsave
from os.path import abspath, dirname
sys.path.append("/groups/multifish/multifish/air_localize_automation/")
import n5_metadata_utils as n5mu
import concurrent.futures

spot_path=sys.argv[1]

# In[3]:


mask=imread("/nrs/multifish/Yuhan/LHA3/segmentation/mask.tif")
print(mask.shape)


# In[2]:


#lb=imread("/Volumes/multifish/Yuhan/LHA3/segmentation/R3_label_affinity.tif")
#roi = np.max(lb)
#print(roi)
#list=np.unique(lb[mask==0])
#print(len(list))
##this is very slow too

#for i in list:
#    lb[lb==i]=0

#lb_bi=np.where(np.isin(lb,list),0,1)  #####this is too slow


# In[4]:


# get appropriate image metadata    
img='/nrs/multifish/Yuhan/LHA3/stitch/R3_LHA3/export.n5/'
subpath='c2/s2'
#voxel size in µm (x, y, z)
vox= n5mu.read_voxel_spacing(img, subpath)
#image size in pixel (x, y, z)
grid=n5mu.read_voxel_grid(img, subpath)
size=grid*vox
print(vox)
print(grid)
print(size)


# In[71]:


b=np.array([-1, -1, -1,-1])
fx=sorted(glob(spot_path))
for f in fx:
    tile=os.path.basename(f)
    name=tile.split('.')[0]
    print(name)
    spot=np.loadtxt(f,delimiter=',')
    m=len(spot)
    spot[:,:3]=spot[:,:3]/vox
    for j in range(0,m):
        if np.all(np.logical_and(spot[j,:3]<=(grid-1),spot[j,:3]>=np.array([(0,0,0)]))):
            if mask[int(spot[j,2]),int(spot[j,1]),int(spot[j,0])]==1:
                b=np.vstack((b,spot[j,:])) 
    b=b[(b != -1).all(axis=1),:]      
    np.savetxt('/nrs/multifish/Yuhan/LHA3/analysis/spot/2_inMask/{}_crop.txt'.format(name), b)
    #np.savetxt('/Volumes/multifish/Yuhan/LHA3/analysis/spot/raw/R2_c0_crop.txt', spot)
    print('raw counts:', m)
    print('counts in mask:', len(b))
    print(b.shape)
