# Code to remove ROIs (1) on edges (2) failed to be imaged in certain FISH round(s) 

import os, sys, zarr
sys.path.append('..')
import numpy as np
import pandas as pd
from glob import glob 
from skimage.io import imread, imsave
from os.path import abspath, dirname
from easi_fish import n5_metadata_utils as n5mu
import warnings
warnings.filterwarnings('ignore')

if __name__ == '__main__':
    fix_dir        = sys.argv[1]  # directory to the fixed image
    reg_dir        = sys.argv[1]  # directory to the registered image
    out_dir        = sys.argv[2]  # where the output files should be saved 
    lb_dir         = sys.argv[3]  # directory to the segmentation mask (tif format accepted here)

# get appropriate image data    
subpath='/c2/s2'
#voxel size in µm (x, y, z) (post-expansion)
vox= n5mu.read_voxel_spacing(fix_dir, subpath)
#image size in pixel (x, y, z)
grid=n5mu.read_voxel_grid(fix_dir, subpath)
#image size in physical space (x, y, z) (post-expansion)
size=grid*vox
print('voxel size is:',vox)
print('image size in pixel unit is:',grid)
print('image size in um unit is:',size)

# get appropriate image data
print("loading images...")
fix = zarr.open(store=zarr.N5Store(fix_dir), mode='r')     
img1 = fix[subpath][:, :, :]

reg = zarr.open(store=zarr.N5Store(reg_dir), mode='r')     
img2 = reg[subpath][:, :, :]

print("all images loaded")

mask=np.full((grid[2], grid[1], grid[0]),1)
mask[img2==0]=0


imsave(out_dir+'/mask.tif',mask)
print("mask generated")
print("mask dimension is:", mask.shape)


lb=imread(lb_dir)
roi = np.max(lb)
print(roi)

# # Get list of ROIs that are fully or partially outside the mask 

bad_roi=np.unique(lb[mask==0])
np.save(out_dir+'/bad_roi_list.npy', bad_roi)

print("# of ROIs rejected:", len(bad_roi))


