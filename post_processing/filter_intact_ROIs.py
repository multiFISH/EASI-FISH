# Code to remove ROIs (1) on edges (2) failed to be imaged in certain FISH round(s) 

import os, sys,z5py
import numpy as np
import pandas as pd
from glob import glob 
from skimage.measure import regionprops
from skimage.io import imread, imsave
from os.path import abspath, dirname
import n5_metadata_utils as n5mu
from scipy.spatial import distance




if __name__ == '__main__':
    img_dir        = sys.argv[1]  # directory to the segmentation mask (tif format accepted here)
    spot_dir       = sys.argv[2]  # directory to folder of detected spot (csv format, usually Airlocalize output)
    out_dir        = sys.argv[3]  # directory where the output file should be saved 
    lb_dir         = sys.argv[4]  # directory to segmentation mask

# get appropriate image data    
subpath='c2/s2'
#voxel size in µm (x, y, z)
vox= n5mu.read_voxel_spacing(img_dir, subpath)
#image size in pixel (x, y, z)
grid=n5mu.read_voxel_grid(img_dir, subpath)
size=grid*vox

# get appropriate image data
im = z5py.File(img_dir, use_zarr_format=False)       
print("loading images")

img2 = im['c0/s0'][:, :, :]
img3 = im['c1/s0'][:, :, :]
img4 = im['c2/s0'][:, :, :]
img5 = im['c3/s0'][:, :, :]
img6 = im['c4/s0'][:, :, :]
img7 = im['c5/s0'][:, :, :]
img8 = im['c6/s0'][:, :, :]
img9 = im['c7/s0'][:, :, :]
img10 = im['c8/s0'][:, :, :]

print("all images loaded")

mask=np.full((grid[2], grid[1], grid[0]),1)
mask[img2==0]=0
mask[img3==0]=0
mask[img4==0]=0
mask[img5==0]=0
mask[img6==0]=0
mask[img7==0]=0
mask[img8==0]=0
mask[img9==0]=0
mask[img10==0]=0

print("mask generated")
print("mask dimension:")
print(mask.shape)


lb=imread(lb_dir)
roi = np.max(lb)
print(roi)

# Optional: 
# to remove the last few z-stacks (z>650) due to high background in the red channel.

mask[650:,:,:]=0


# Get list of ROIs that are fully or partially outside the mask 

list=np.unique(lb[mask==0])
print("# of cells rejected:", len(list))

# ####Get metadata for ROI: ID, centroid position, size, distance to (0,0,0) and aspect ratio####

df = pd.DataFrame(np.empty([roi, 0]))
lb_stat = regionprops(lb)
for i in range(0,roi): 
    df.loc[df.index[i], 'roi'] = i+1
    df.loc[df.index[i], 'z'] = lb_stat[i].centroid[0]
    df.loc[df.index[i], 'y'] = lb_stat[i].centroid[1]
    df.loc[df.index[i], 'x'] = lb_stat[i].centroid[2]
    df.loc[df.index[i], 'area'] = lb_stat[i].area
    df.loc[df.index[i], 'Distance'] = distance.euclidean(lb_stat[i].centroid,[0,0,0])
    df.loc[df.index[i], 'minor_axis_length'] = lb_stat[i].minor_axis_length
    df.loc[df.index[i], 'major_axis_length'] = lb_stat[i].major_axis_length
    df.loc[df.index[i], 'aspect_ratio'] = lb_stat[i].minor_axis_length/lb_stat[i].major_axis_length

####Filter out ROIs that 1) not in mask; 2) have high background in channel 546

df_filtered=df.loc[~df['roi'].isin(list)]

df_filtered.to_csv(out_dir, index=False)





