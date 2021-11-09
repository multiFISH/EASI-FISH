
import os
import sys
import numpy as np
import pandas as pd
from glob import glob
from os.path import abspath, dirname
from skimage.measure import regionprops

def roi_prop(lb,s=[0.92,0.92,0.84],ex=2):

    """
    Returns ROI position, volume, aspect ratio.  
    This function uses the regionprops function implemented in skimage  
    lb: segmentation mask
    s: pixel size in µm for [x,y,z],default [0.92,0.92,0.84], 
    ex: linear expansion factor, default is 2
    """
    lb_id = np.unique(lb[lb != 0])
    roi=len(lb_id)
    df = pd.DataFrame(np.empty([roi, 0]))
    lb_stat = regionprops(lb)

    for i in range(0,roi): 
        df.loc[df.index[i], 'roi'] = i+1
        df.loc[df.index[i], 'z'] = lb_stat[i].centroid[0]*s[2]/ex
        df.loc[df.index[i], 'y'] = lb_stat[i].centroid[1]*s[1]/ex
        df.loc[df.index[i], 'x'] = lb_stat[i].centroid[2]*s[0]/ex
        df.loc[df.index[i], 'area'] = lb_stat[i].area*s[2]*s[1]*s[0]/(ex**3)
        df.loc[df.index[i], 'minor_axis_length'] = lb_stat[i].minor_axis_length
        df.loc[df.index[i], 'major_axis_length'] = lb_stat[i].major_axis_length
        df.loc[df.index[i], 'eccentricity'] = lb_stat[i].minor_axis_length/lb_stat[i].major_axis_length
        df.loc[df.index[i], 'solidity'] = lb_stat[i].solidity
        
    return(df)

