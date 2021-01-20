#!/usr/bin/env python
# coding: utf-8
# ####New lipofuscin subtraction method test###

import os, sys
import numpy as np
import pandas as pd
from glob import glob 
from skimage.measure import regionprops
from skimage.io import imread, imsave
from os.path import abspath, dirname
from scipy.spatial import distance

if __name__ == '__main__':
    
    spotcount_dir  = sys.argv[1]  # directory to spot count/neuron (csv format)
    lipocount_dir  = sys.argv[2]  # directory to estimated lipofuscin count/neuron (csv format)
    roi_dir        = sys.argv[4]  # directory to roi metadata
    output_dir     = sys.argv[5]  # directory where output should be stored


spots_raw=pd.read_csv(spotcount_dir,index_col=0, sep=',')
lipo=pd.read_csv(lipocount_dir,index_col=0, sep=',')
roi=pd.read_csv(roi_dir,sep=',', index_col=0)

spots_raw=spots_raw[spots_raw.index.isin(roi.index)]
lipo=lipo[lipo.index.isin(roi.index)]

lipo['c0_mean']=lipo[['R4_c0_lipo','R5_c0_lipo','R6_c0_lipo','R7_c0_lipo','R8_c0_lipo','R9_c0_lipo','R10_c0_lipo']].mean(axis=1)
lipo['c0_median']=lipo[['R4_c0_lipo','R5_c0_lipo','R6_c0_lipo','R7_c0_lipo','R8_c0_lipo','R9_c0_lipo','R10_c0_lipo']].median(axis=1)

lipo['c1_mean']=lipo[['R4_c0_lipo','R5_c0_lipo','R6_c0_lipo','R7_c0_lipo','R8_c0_lipo','R9_c0_lipo','R10_c0_lipo']].mean(axis=1)
lipo['c1_median']=lipo[['R4_c0_lipo','R5_c0_lipo','R6_c0_lipo','R7_c0_lipo','R8_c0_lipo','R9_c0_lipo','R10_c0_lipo']].median(axis=1)

# For R2 and R3, subtract lipofuscin spots in cells based on the identified lipofuscins in corresponding rounds; 
# For R4-R10, if the raw spots > 500 in either c0 or c1, then use the median lipofuscin counts (to avoid subtracting real spots); if the raw spots <500, then use the identified lipofuscins in corresponding rounds;

# Handling of R1, R2, R3

spots_new=spots_raw.copy()
for i in roi.index:
    spots_new.loc[i,'R1_c0']=spots_new.loc[i,'R1_c0']-lipo.loc[i,'R1_c0_lipo']
    spots_new.loc[i,'R1_c1']=spots_new.loc[i,'R1_c1']-lipo.loc[i,'R1_c1_lipo']
    spots_new.loc[i,'R2_c0']=spots_new.loc[i,'R2_c0']-lipo.loc[i,'R2_c0_lipo']
    spots_new.loc[i,'R2_c1']=spots_new.loc[i,'R2_c1']-lipo.loc[i,'R2_c1_lipo']
    spots_new.loc[i,'R3_c0']=spots_new.loc[i,'R3_c0']-lipo.loc[i,'R3_c0_lipo']
    spots_new.loc[i,'R3_c1']=spots_new.loc[i,'R3_c1']-lipo.loc[i,'R3_c1_lipo']

# Handling of R4-R10
spots_new=spots_raw.copy()
Uplimit=200
for i in roi.index:
    if (spots_new.loc[i,'R4_c0']>Uplimit or spots_new.loc[i,'R4_c1']>Uplimit)==True:
        spots_new.loc[i,'R4_c0']=spots_new.loc[i,'R4_c0']-lipo.loc[i,'c0_median']
        spots_new.loc[i,'R4_c1']=spots_new.loc[i,'R4_c1']-lipo.loc[i,'c1_median']
    else:
        spots_new.loc[i,'R4_c0']=spots_new.loc[i,'R4_c0']-lipo.loc[i,'R4_c0_lipo']
        spots_new.loc[i,'R4_c1']=spots_new.loc[i,'R4_c1']-lipo.loc[i,'R4_c1_lipo']

for i in roi.index:
    if (spots_new.loc[i,'R5_c0']>Uplimit or spots_new.loc[i,'R5_c1']>Uplimit)==True:
        spots_new.loc[i,'R5_c0']=spots_new.loc[i,'R5_c0']-lipo.loc[i,'c0_median']
        spots_new.loc[i,'R5_c1']=spots_new.loc[i,'R5_c1']-lipo.loc[i,'c1_median']
    else:
        spots_new.loc[i,'R5_c0']=spots_new.loc[i,'R5_c0']-lipo.loc[i,'R5_c0_lipo']
        spots_new.loc[i,'R5_c1']=spots_new.loc[i,'R5_c1']-lipo.loc[i,'R4_c1_lipo']

for i in roi.index:
    if (spots_new.loc[i,'R6_c0']>Uplimit or spots_new.loc[i,'R6_c1']>Uplimit)==True:
        spots_new.loc[i,'R6_c0']=spots_new.loc[i,'R6_c0']-lipo.loc[i,'c0_median']
        spots_new.loc[i,'R6_c1']=spots_new.loc[i,'R6_c1']-lipo.loc[i,'c1_median']
    else:
        spots_new.loc[i,'R6_c0']=spots_new.loc[i,'R6_c0']-lipo.loc[i,'R6_c0_lipo']
        spots_new.loc[i,'R6_c1']=spots_new.loc[i,'R6_c1']-lipo.loc[i,'R6_c1_lipo']

for i in roi.index:
    if (spots_new.loc[i,'R7_c0']>Uplimit or spots_new.loc[i,'R7_c1']>Uplimit)==True:
        spots_new.loc[i,'R7_c0']=spots_new.loc[i,'R7_c0']-lipo.loc[i,'c0_median']
        spots_new.loc[i,'R7_c1']=spots_new.loc[i,'R7_c1']-lipo.loc[i,'c1_median']
    else:
        spots_new.loc[i,'R7_c0']=spots_new.loc[i,'R7_c0']-lipo.loc[i,'R7_c0_lipo']
        spots_new.loc[i,'R7_c1']=spots_new.loc[i,'R7_c1']-lipo.loc[i,'R7_c1_lipo']

for i in roi.index:
    if (spots_new.loc[i,'R8_c0']>Uplimit or spots_new.loc[i,'R8_c1']>Uplimit)==True:
        spots_new.loc[i,'R8_c0']=spots_new.loc[i,'R8_c0']-lipo.loc[i,'c0_median']
        spots_new.loc[i,'R8_c1']=spots_new.loc[i,'R8_c1']-lipo.loc[i,'c1_median']
    else:
        spots_new.loc[i,'R8_c0']=spots_new.loc[i,'R8_c0']-lipo.loc[i,'R8_c0_lipo']
        spots_new.loc[i,'R8_c1']=spots_new.loc[i,'R8_c1']-lipo.loc[i,'R8_c1_lipo']

for i in roi.index:
    if (spots_new.loc[i,'R9_c0']>Uplimit or spots_new.loc[i,'R9_c1']>Uplimit)==True:
        spots_new.loc[i,'R9_c0']=spots_new.loc[i,'R9_c0']-lipo.loc[i,'c0_median']
        spots_new.loc[i,'R9_c1']=spots_new.loc[i,'R9_c1']-lipo.loc[i,'c1_median']
    else:
        spots_new.loc[i,'R9_c0']=spots_new.loc[i,'R9_c0']-lipo.loc[i,'R9_c0_lipo']
        spots_new.loc[i,'R9_c1']=spots_new.loc[i,'R9_c1']-lipo.loc[i,'R9_c1_lipo']

for i in roi.index:
    if (spots_new.loc[i,'R10_c0']>Uplimit or spots_new.loc[i,'R10_c1']>Uplimit)==True:
        spots_new.loc[i,'R10_c0']=spots_new.loc[i,'R10_c0']-lipo.loc[i,'c0_median']
        spots_new.loc[i,'R10_c1']=spots_new.loc[i,'R10_c1']-lipo.loc[i,'c1_median']
    else:
        spots_new.loc[i,'R10_c0']=spots_new.loc[i,'R10_c0']-lipo.loc[i,'R10_c0_lipo']
        spots_new.loc[i,'R10_c1']=spots_new.loc[i,'R10_c1']-lipo.loc[i,'R10_c1_lipo']


spots_new.to_csv(output_dir)





