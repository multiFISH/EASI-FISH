# Spot counts for cells with highly expressed genes (dense spots)
# 1. Measure total intensity of every ROI after bleed-through correction and background subtraction.
# 2. Calculate the number of spot from total intensity based on unitary spot fluorescence intensity
# 3. Correlate the number of spots (from air-localize) with the total fluorescence intensity/voxel in each ROI and determine a 'cutoff'. 
#    Spot count > cutoff: use spot count converted based on total fluorescence intensity; 
#    Spot count < cutoff: use spot count from Airlocalize

import os, sys,z5py
import numpy as np
import pandas as pd
from glob import glob 
from skimage.measure import regionprops
from skimage.io import imread, imsave
from os.path import abspath, dirname
from scipy import stats
from scipy.stats import skewnorm,lognorm
from scipy.optimize import minimize



if __name__ == '__main__':
    spotcount_dir      = sys.argv[1]  # directory to assigned spots per neuron based on airlocalize (csv format)
    roi_dir            = sys.argv[2]  # directory to file containing the ROI metadata (neuron volume, etc.)
    GeneName_dir       = sys.argv[3]  # directory to file containing marker gene metadata (gene name, corresponding imaging round and image channel)
    spot_dir           = sys.argv[4]  # directory to folder of airlocalize output (1 file/gene, txt format)
    intensity_dir      = sys.argv[5]  # directory to folder of intensity measurement output (1 file/gene, csv format)
    output_dir         = sys.argv[6]  # directory where output should be stored

spotcount=pd.read_csv(spotcount_dir,sep=',', index_col=0)
roi=pd.read_csv(roi_dir,sep=',', index_col=0)
GeneName=pd.read_csv(GeneName_dir,sep=',', index_col=0)

### Identify unitary spot fluorescence intensity for every gene
fx=sorted(glob(spot_dir+'/*.txt'))
for f in fx:
    r=os.path.basename(f).split('.')[0]
    spot=np.loadtxt(f, delimiter=',')
    vox=[0.92,0.92,0.84]
    spot[:,:3]=spot[:,:3]/vox  # convert from physical unit to pixel unit
    for i in range(2):
        spot=spot[np.logical_and(spot[:,i]<=1500,spot[:,i]>250)]
    spot=spot[np.logical_and(spot[:,2]<=650,spot[:,2]>150)]   ##remove spots on edges (eliminate false detection)
    spot_int= spot[:,3]
    spot_int=spot_int[spot_int!=-8.0]
    n,b=np.histogram(spot_int, bins=5000)
    GeneName.loc['%s' % (r), 'single_spot_intensity']=b[np.argwhere(n == n.max())][0][0] 
    ##Note that the histogram maximum is used as an estimate for single spot intensity. We also tried fitting the data to a skewed normal (or log-normal) distribution and then estimate the peak (see below).  
#     ae, loce, scalee = skewnorm.fit(spot_int)
#     def skew_fit(n):
#         return -skewnorm.pdf(n, ae, loce, scalee)
#     GeneName.loc['%s' % (r), 'single_spot_intensity']=minimize(skew_fit,0,method='Powell').x

# df_mean is mean_fluorescence_intensity (after background subtraction) 
# df_total is total_fluorescence_intensity (after subtracting background)
# df_count is spot count calculated from total fluorescence intensity
df_mean = pd.DataFrame(data=np.empty([len(roi),0]), index=roi.index, dtype=float)
df_total = pd.DataFrame(data=np.empty([len(roi),0]), index=roi.index, dtype=float)
df_count = pd.DataFrame(data=np.empty([len(roi),0]), index=roi.index, dtype=float)
fx=sorted(glob(intensity_dir+"/*_c0_intensity.csv"))
for f in fx:
    r=os.path.basename(f).split('_')[0]
    c=os.path.basename(f).split('_')[1]
    cell_int=pd.read_csv(f,sep=',', index_col=0)
    cell_int=cell_int[cell_int.index.isin(roi.index)]  ## only include intact ROIs###
    n,b=np.histogram(cell_int['mean_intensity'], bins=1000) ## Idenfity background###
    bg=b[np.argwhere(n == n.max())][0][0]                   ## Idenfity background###
    df_mean['%s_%s' % (r,c)]=np.maximum(0,cell_int['mean_intensity']-bg)
    df_total['%s_%s' % (r,c)]=np.maximum(0,(cell_int['mean_intensity']-bg))*roi['area']
    df_count['%s_%s' % (r,c)]=df_total['%s_%s' % (r,c)]/GeneName.loc['%s_%s' % (r,c), 'single_spot_intensity']

fx=sorted(glob(intensity_dir+"/*_c1_intensity.csv"))
for f in fx:
    r=os.path.basename(f).split('_')[0]
    c=os.path.basename(f).split('_')[1]
    cell_int=pd.read_csv(f,sep=',', index_col=0)
    cell_int=cell_int[cell_int.index.isin(roi.index)]  ## only include intact ROIs###
    n,b=np.histogram(cell_int['mean_intensity'], bins=1000) ## Idenfity background###
    bg=b[np.argwhere(n == n.max())][0][0]                   ## Idenfity background###
    df_mean['%s_%s' % (r,c)]=np.maximum(0,cell_int['mean_intensity']-bg)
    df_total['%s_%s' % (r,c)]=np.maximum(0,(cell_int['mean_intensity']-bg))*roi['area']
    df_count['%s_%s' % (r,c)]=df_total['%s_%s' % (r,c)]/GeneName.loc['%s_%s' % (r,c), 'single_spot_intensity']

fx=sorted(glob(intensity_dir+"/*_c3_intensity.csv"))
for f in fx:
    r=os.path.basename(f).split('_')[0]
    c=os.path.basename(f).split('_')[1]
    cell_int=pd.read_csv(f,sep=',', index_col=0)
    cell_int=cell_int[cell_int.index.isin(roi.index)]   ## only include intact ROIs###
    df_mean['%s_%s' % (r,c)]=np.maximum(0,cell_int['mean_intensity'])   #No background subtraction for c3, because it has been handled in the intensity measurement step
    df_total['%s_%s' % (r,c)]=np.maximum(0,(cell_int['mean_intensity']))*roi['area']
    df_count['%s_%s' % (r,c)]=df_total['%s_%s' % (r,c)]/GeneName.loc['%s_%s' % (r,c), 'single_spot_intensity']

df_cutoff = spotcount.copy()
for i in df_count.columns:
    density=spotcount[i]/(roi['area']*2*2*2/(0.92*0.92*0.84) # convert um^3 to voxel values
    for j in density[density>0.01].index:  ##this threshold corresponds to spot-spot distance ~1.3 um apart
        df_cutoff.loc[j,i]=df_count.loc[j,i].copy()
        
df_cutoff.to_csv(output_dir+'/spotcount_dense_spot_corrected.csv')
