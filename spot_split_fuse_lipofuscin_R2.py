import numpy as np
from skimage.io import imread, imsave
from glob import glob 
import os, sys
import pandas as pd
from skimage.io import imread, imsave
from os.path import abspath, dirname
sys.path.append("/groups/multifish/multifish/air_localize_automation/")
import n5_metadata_utils as n5mu
from scipy.spatial import distance

if __name__ == "__main__":

    spot_c0_path      = sys.argv[1]
    spot_c1_path      = sys.argv[2]
    outdir            = sys.argv[3]
    round             = sys.argv[4]
    
# get appropriate image data    
img='/nrs/multifish/Yuhan/LHA3/stitch/R3_LHA3/export.n5/'
subpath='c2/s0'
subpath_2='c2/s2'
#voxel size in µm (x, y, z)
vox= n5mu.read_voxel_spacing(img, subpath)
vox2= n5mu.read_voxel_spacing(img, subpath_2)
scale=vox2/vox
#image size in pixel (x, y, z)
grid=n5mu.read_voxel_grid(img, subpath)
shape=grid*vox
print(vox)
print(grid)


spot_c0=np.loadtxt(spot_c0_path)
spot_c1=np.loadtxt(spot_c1_path,delimiter=',')
print('total number of spots in c0')
print(spot_c0.shape)
print('total number of spots in c1')
print(spot_c1.shape)


def spot_trim(spot,d,a,b):
    spot=spot[np.logical_and(spot[:,d]<=b,spot[:,d]>a)]
    return spot

def start_points(size, split_size, overlap):
    points = [0]
    stride = int(split_size * (1-overlap))
    counter = 1
    while True:
        pt = stride * counter
        if pt + split_size >= size:
            points.append(size-split_size)
            break
        else:
            points.append(pt)
        counter += 1
    points=np.ceil(points).astype(int)
    return points

def init_points(corner, size, split_size, overlap):
        edge=split_size*(overlap*0.5)
        if corner==0:
            start=corner
            end=corner+split_size-edge
        else:
            if corner+split_size>=size:
                rest=(size-overlap*split_size)%(split_size*(1-overlap))
                start=size-rest-edge
                end=size
            else: 
                start=corner+edge
                end=corner+split_size-edge
        return start, end

split_x = 1000
split_y = 1000
split_z = 200
X_points = start_points(grid[0], split_x, 0.1)
Y_points = start_points(grid[1], split_y, 0.1)
Z_points = start_points(grid[2], split_z, 0.1)
print(X_points,Y_points,Z_points)


#split image and save
print('splitting---')
d1={}
d2={}
d3={}
d4={}

spot_c0[:,:3]=spot_c0[:,:3]/vox
spot_c1[:,:3]=spot_c1[:,:3]/vox

for i in X_points :
     for j in Y_points :
        for k in Z_points :
            c0=spot_trim(spot_c0,2,k,k+split_z)
            c1=spot_trim(spot_c1,2,k,k+split_z)
            c0=spot_trim(spot_c0,1,j,j+split_y)
            c1=spot_trim(spot_c1,1,j,j+split_y)
            c0=spot_trim(spot_c0,0,i,i+split_x)
            c1=spot_trim(spot_c1,0,i,i+split_x)
            dist=distance.cdist(c0[:,:3], c1[:,:3])
            lipo_c0=np.argwhere(dist<=3)[:,0]
            lipo_c1=np.argwhere(dist<=3)[:,1]
            d1['lipo_%d_%d_%d' % (i,j,k)]=c0[lipo_c0,:]
            d2['lipo_%d_%d_%d' % (i,j,k)]=c1[lipo_c1,:]
            d3['real_%d_%d_%d' % (i,j,k)]=np.delete(c0, lipo_c0, axis=0)
            d4['real_%d_%d_%d' % (i,j,k)]=np.delete(c1, lipo_c1, axis=0)
            print(i,j,k, c0.shape)

# In[61]:


print('merging lipofuscin spots')
a=[-1,-1,-1,-1]
for d in d1:
    i=int(d.split('_')[1])
    j=int(d.split('_')[2])
    k=int(d.split('_')[3])
    spot=d1[d]
    z1,zN=init_points(k, grid[2], split_z, overlap=0.1)
    y1,yN=init_points(j, grid[1], split_y, overlap=0.1)
    x1,xN=init_points(i, grid[0], split_x, overlap=0.1)
    spot=spot[np.logical_and(spot[:,0]<=xN,spot[:,0]>x1)]
    spot=spot[np.logical_and(spot[:,1]<=yN,spot[:,1]>y1)]
    spot=spot[np.logical_and(spot[:,2]<=zN,spot[:,2]>z1)]
    a=np.vstack((a,spot))
    a=a[(a != -1).all(axis=1),:]  
    
a[:,:3]=a[:,:3]/scale    
np.savetxt(outdir+'/{}_c0_lipo.txt'.format(round), a)
print('lipo_c0 saved')
print(a.shape)

# In[47]:


b=[-1,-1,-1,-1]
for d in d2:
    i=int(d.split('_')[1])
    j=int(d.split('_')[2])
    k=int(d.split('_')[3])
    spot=d2[d]
    z1,zN=init_points(k, grid[2], split_z, overlap=0.1)
    y1,yN=init_points(j, grid[1], split_y, overlap=0.1)
    x1,xN=init_points(i, grid[0], split_x, overlap=0.1)
    spot=spot[np.logical_and(spot[:,0]<=xN,spot[:,0]>x1)]
    spot=spot[np.logical_and(spot[:,1]<=yN,spot[:,1]>y1)]
    spot=spot[np.logical_and(spot[:,2]<=zN,spot[:,2]>z1)]
    b=np.vstack((b,spot))

b=b[(b != -1).all(axis=1),:]  
b[:,:3]=b[:,:3]/scale
np.savetxt(outdir+'/{}_c1_lipo.txt'.format(round), b)    
print('lipo_c1 saved')
print(b.shape)

# In[53]:


print('subtracting lipofuscin spots')
c=[-1,-1,-1,-1]
for d in d3:
    i=int(d.split('_')[1])
    j=int(d.split('_')[2])
    k=int(d.split('_')[3])
    spot=d3[d]
    z1,zN=init_points(k, grid[2], split_z, overlap=0.1)
    y1,yN=init_points(j, grid[1], split_y, overlap=0.1)
    x1,xN=init_points(i, grid[0], split_x, overlap=0.1)
    spot=spot[np.logical_and(spot[:,0]<=xN,spot[:,0]>=x1)]
    spot=spot[np.logical_and(spot[:,1]<=yN,spot[:,1]>=y1)]
    spot=spot[np.logical_and(spot[:,2]<=zN,spot[:,2]>=z1)]
    c=np.vstack((c,spot))

c=c[(c != -1).all(axis=1),:]  
c[:,:3]=c[:,:3]/scale  
np.savetxt(outdir+'/{}_c0_true.txt'.format(round), c)
print('real_c0 saved')
print(c.shape) 

# In[50]:


e=[-1,-1,-1,-1]
for d in d4:
    i=int(d.split('_')[1])
    j=int(d.split('_')[2])
    k=int(d.split('_')[3])
    spot=d4[d]
    z1,zN=init_points(k, grid[2], split_z, overlap=0.1)
    y1,yN=init_points(j, grid[1], split_y, overlap=0.1)
    x1,xN=init_points(i, grid[0], split_x, overlap=0.1)
    spot=spot[np.logical_and(spot[:,0]<=xN,spot[:,0]>=x1)]
    spot=spot[np.logical_and(spot[:,1]<=yN,spot[:,1]>=y1)]
    spot=spot[np.logical_and(spot[:,2]<=zN,spot[:,2]>=z1)]
    e=np.vstack((e,spot))


e=e[(e != -1).all(axis=1),:]  
e[:,:3]=e[:,:3]/scale      
np.savetxt(outdir+'/{}_c1_true.txt'.format(round), e)
print('real_c1 saved')
print(e.shape)

print(spot_c0.shape)
print(spot_c1.shape)


