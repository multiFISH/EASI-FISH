# Code used to flag over-segmentation errors
# A. maximize the number of identified errors
# B. minimize false detection of well-segemented cells 
# 
# The over-segmentation pairs are identified with the following criteria: 
# 1) the two segments have a very high correlation (corr > 0.99) 
# 
# 2) the two segments have to be less than 55 pixels apart(centroid position distance)
# 
# 3) at least one of the segments have to be bigger than 7000 pixels in size
# 
# 4) The two segments have to be touching each other
# 
# 5) if both segments are bigger than 20000 in size, it will be flagged and manually inspected
# 
# 6) for ROIs that have more than two corresponding pairs, the corresponding segments will be ranked by correlation coefficient, and they will be flagged for manual inspection
# 
# #1 and #2 are for initial identification of oversegmentation errors 
# 
# #3, #4, #5 are for eliminating the false detections  

import os, sys
import numpy as np
import pandas as pd
from glob import glob 
from skimage.io import imread, imsave
from os.path import abspath, dirname
from scipy import stats
from scipy.spatial import distance,cKDTree
from skimage.measure import regionprops

if __name__ == '__main__':
    count_dir      = sys.argv[1]  # directory to spot count/neuron (csv format)
    roi_dir        = sys.argv[2]  # directory to neuron metadata (csv format)
    lb_dir         = sys.argv[3]  # directory to segmentation mask (tif format)
    out_dir        = sys.argv[4]  # directory to save output

df=pd.read_csv(count_dir,sep=',', index_col=0)
roi=pd.read_csv(roi_dir,sep=',',index_col=0)

###Get correlation matrix
corr_raw =df.T.corr()
s_raw = corr_raw.stack()
ii_raw = s_raw[np.logical_and(s_raw > 0.998, s_raw<1.0)].index.tolist()
test=np.asarray(ii_raw)
test.sort(axis=1)
test=np.unique(test,axis=0)

cand={}
for i in range(0,len(test)):
    a=roi.loc[float(test[i,0])].to_numpy()[:3]
    b=roi.loc[float(test[i,1])].to_numpy()[:3]
    dist=distance.euclidean(a,b)
    if dist<55 and dist>0:
        a_area=roi.loc[float(test[i,0])]['area']
        b_area=roi.loc[float(test[i,1])]['area']
        if np.maximum(a_area,b_area) >7000:
            c=corr_raw.loc[test[i,0],test[i,1]]
            cand[i] = np.append(test[i,:],c)
            cand[i] = np.append(cand[i],dist)
            if np.minimum(a_area,b_area) >20000:
                cand[i]=np.append(cand[i], str('atten'))
            else:
                cand[i]=np.append(cand[i], str('--'))

m=pd.DataFrame.from_dict(data=cand, orient='index')
m = m.rename(columns={1:'cell_A', 2:'cell_B', 3: 'corr', 4:'dist',5: 'min_size_20000'})
m['cell_A'], m['cell_B'] = np.where(m['cell_A'] > m['cell_B'], [m['cell_B'],m['cell_A'] ], [m['cell_A'] , m['cell_B']])

lb=imread(lb_dir)
lb_stat=regionprops(lb)

select={}
for k in range(0,len(m)):
    #a=np.argwhere(lb_4==n.iloc[k]['cell_A'])
    #b=np.argwhere(lb_4==n.iloc[k]['cell_B'])
    id1=m.iloc[k]['cell_A'].copy()
    id2=m.iloc[k]['cell_B'].copy()
    a=lb_stat[int(id1-1)].coords
    b=lb_stat[int(id2-1)].coords
    kdtree_a = cKDTree(a)
    kdtree_b = cKDTree(b)
    #neighbors = kdtree_a.query_ball_tree(kdtree_b, 1)
    nnn=kdtree_a.count_neighbors(kdtree_b,1)
    if nnn>0:
        select[k]=m.iloc[k]
        if nnn<50:
            select[k]=np.append(select[k], str('less_50'))
        else:
            select[k]=np.append(select[k], str('--'))
select=pd.DataFrame.from_dict(data=select, orient='index')
select = select.rename(columns={0:'cell_A', 1:'cell_B', 2: 'corr', 3:'dist',4: 'min_size_20000',5:'touch'})

select.to_csv(out_dir+/'flag_oversegmentation.csv')
## This output a list of oversegmented ROI pairs, with some flagged for manual inspection. After manual inspection, the list can be used to merge the oversegmented ROI pairs.   
