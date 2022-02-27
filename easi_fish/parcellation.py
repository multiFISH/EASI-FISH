import os, sys
import numpy as np
import pandas as pd
from glob import glob 
from scipy import spatial
import matplotlib.pyplot as plt
from matplotlib.colors import DivergingNorm
import math

def relative_expression(roi, spotcount, gene1, gene2, radius): 
    """
    This is to compute the relative spatial expression relationship of two genes, gene1 and gene2. 
    
    roi: pandas dataframe containing neuron id and their x,y,z positions (in the first 3 columns)
    spotcount: pandas dataframe containing neuron id and their gene expression (spot count)
    gene1 and gene2: genes to compute the relative relationship
    radius: the neighborhood radius used to compute the spatial gene expression 

    """   
    gene1=str(gene1)
    gene2=str(gene2)
    roi.loc[:, 'cluster']='N'
    roi.loc[spotcount[spotcount[gene1]>50].index, 'cluster']=gene1
    roi.loc[spotcount[spotcount[gene2]>50].index, 'cluster']=gene2
    roi.loc[spotcount[(spotcount[gene1]>50)&(spotcount[gene2]>50)].index, 'cluster']=gene1+'_'+gene2
    X=roi.to_numpy()[:,:3] 
    neuron=spatial.KDTree(X)
    neighbors=neuron.query_ball_point(X,radius)
    roi['cluster']=roi['cluster'].astype(str)
    roi['fraction_%s_%s' % (str(gene1), str(gene2))]=0

    ind1=roi.columns.get_loc('cluster')
    ind2=roi.columns.get_loc('fraction_%s_%s' % (str(gene1), str(gene2)))
    for i in range(0,len(neighbors)):
        x=[]
        for j in neighbors[i]:
            x=np.append(x, roi.iloc[j,ind1])
        a,b=np.unique(x,return_counts=True)
        if np.any(a==gene1) and np.any(a==gene2):
            c=b[np.argwhere(a==gene1)]
            d=b[np.argwhere(a==gene2)]
            roi.iloc[i,ind2]=float((c-d)/(c+d))
        else:
            if np.any(a==gene1):
                roi.iloc[i,ind2]=1
            if np.any(a==gene2):
                roi.iloc[i,ind2]=-1
    return roi

def spatial_enrichment(roi, spotcount, gene, radius): 
    """
    This is to compute the spatial enrichment of a gene. 
    
    roi: pandas dataframe containing neuron id and their x,y,z positions (in the first 3 columns)
    spotcount: pandas dataframe containing neuron id and their gene expression (spot count)
    gene: gene to compute the spatial expression
    radius: the neighborhood radius used to compute the spatial gene expression 

    """   
	gene=str(gene)
	roi.loc[:, 'cluster']='N'
	roi.loc[spotcount[spotcount[gene1]>50].index, 'cluster']=gene

	X=roi.to_numpy()[:,:3] 
	neuron=spatial.KDTree(X)
	neighbors=neuron.query_ball_point(X,radius)
	roi['cluster']=roi['cluster'].astype(str)
	roi['fraction_%s' % (str(gene))]=0

	ind1=roi.columns.get_loc('cluster')
	ind2=roi.columns.get_loc('fraction_%s' % (str(gene)))
	for i in range(0,len(neighbors)):
		x=[]
		for j in neighbors[i]:
			x=np.append(x, roi.iloc[j,ind1])
		a,b=np.unique(x,return_counts=True)
		if np.any(a==gene1):
			c=b[np.argwhere(a==gene)]
			roi.iloc[i,ind2]=float(c/len(x))
    return roi

def plot_relative_expression(roi, column, num_z, invert_x=False, invert_y=False, invert_z=False):
    """
    Plot the relative spatial expression relationship of two genes, as computed with the 'relative_expression' function. 
    
    roi: pandas dataframe containing neuron id and their x,y,z positions (in the first 3 columns)
    column: string, the column name containing the relative expression 
    num_z: number of axial levels to plot.   
    invert_x, invert_y, invert_z: True or False, whether to invert the x,y and z axis. Default is False
    
    """   
    width_ratio=np.ones(num_z)
    width_ratio=np.append(width_ratio,0.1)               
    fig,ax=plt.subplots(1,num_z+1,figsize=(num_z*5+1,5),dpi=150,gridspec_kw={"width_ratios":width_ratio})
    ind=0
    A=roi.copy()
    s=round(math.ceil(roi.z.max())/num_z,-1)
    column=str(column)
    if invert_z:
        A.z=A.z.max()-A.z
    for n in range(1,num_z+1):
        B=A[(A.z>((n-1)*s))&(A.z<=(n*s))]
        x=roi[roi.index.isin(B.index)][column].astype(float)
        a=ax.flatten()[ind].scatter(B.to_numpy()[:,2],(B.to_numpy()[:,1]), c=x, 
                                    norm=DivergingNorm(vcenter=0),
                                    marker='o', s=10, cmap=plt.cm.coolwarm,alpha=1)
        ax.flatten()[ind].yaxis.set_tick_params(pad=10,labelsize=12)
        ax.flatten()[ind].xaxis.set_tick_params(pad=3,labelsize=12)
        ax.flatten()[ind].set_xlim(0,850)
        ax.flatten()[ind].set_ylim(0,850)
        ax.flatten()[ind].set(adjustable='box', aspect='equal')
        ax.flatten()[ind].set_title(str(s*(n-1))+'-'+str(s*n)+'µm (z axis)',fontsize=16)
        ax.flatten()[ind].xaxis.set_ticklabels([])
        ax.flatten()[ind].yaxis.set_ticklabels([])
        ax.flatten()[ind].xaxis.set_ticks([])
        ax.flatten()[ind].yaxis.set_ticks([])
        if invert_y: 
            ax.flatten()[ind].invert_yaxis()
        if invert_x:
            ax.flatten()[ind].invert_xaxis()
        ind = ind + 1
    cb=fig.colorbar(a, cax=ax.flatten()[ind], shrink=0.1,aspect=10,pad=0.5)
    cb.ax.set_title(column,size=12)
    plt.tight_layout()






