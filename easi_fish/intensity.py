# Code to measure mean fluorescence intensity (instead of spot count) from each ROI. 
# The integrated intensity = mean intensity * size of ROI

import os, sys, zarr
import numpy as np
import pandas as pd
from glob import glob
from skimage.measure import regionprops
from skimage.io import imread, imsave

def measure_intensity(lb, img_dir, channel, scale='s2'):
	"""
    Returns spot counts for each ROI. 
    
    lb: segmenntation mask
    img_dir: n5 image file to measure fluorescence intensity
	channel: image channnel, e.g. c0, c1
	scale: image channel, e.g. s0, s1. Default to 's2', 4x4x2 downsampled 
    """
	lb_id = np.unique(lb[lb!= 0])
	roi=len(lb_id)
	if type(img_dir)==str:
		# get n5 image data
		im = zarr.open(store=zarr.N5Store(img_dir), mode='r')     
		img = im[channel+'/'+scale][:,:,:] 
		if channel == 'c3': #default c3 to the channel where bleed-through correction is needed, modify if it does not apply. 
			dapi=im['c2/'+scale][...] #default c2 to DAPI channel, change if dapi is in a different channel
			lo=np.percentile(np.ndarray.flatten(dapi),99.5)
			bg_dapi=np.percentile(np.ndarray.flatten(dapi[dapi!=0]),1)
			bg_img=np.percentile(np.ndarray.flatten(img[img!=0]),1)
			dapi_factor=np.median((img[dapi>lo] - bg_img)/(dapi[dapi>lo] - bg_dapi))
			img = np.maximum(0, img - bg_img - dapi_factor * (dapi - bg_dapi)).astype('float32')
			print('bleed_through:',dapi_factor)
			print('DAPI background:',bg_dapi)
			print('c3 background:',bg_img)
		df_mean = pd.DataFrame(data=np.empty([roi,1]), columns=['roi'], dtype=object)
		lb_stat = regionprops(lb,intensity_image=img)
		for i in range(0,roi):
			df_mean.loc[i, 'roi'] = lb_stat[i].label
			df_mean.loc[i, channel] = lb_stat[i].mean_intensity
	else:
		print('Image file is in wrong format')
	
	return(df_mean)

		
