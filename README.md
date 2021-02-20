# EASI-FISH analysis toolbox
Table of Contents
=================
   * [Description](#description)
   * [Pipeline](#pipeline)
   * [Modules](#modules)
      * [Stitching](#stitching)
      * [Registration](#registration)
      * [Segmentation](#segmentation)
      * [Spot detection](#spot-detection)
   * [Additional information](#additional-information)
      * [Post processing](#post-processing)  
      * [Example data](#example-data)

## Description #
This pipeline is used to analyze large-scale, multi-round, high-resolution image data acquired using EASI-FISH (Expansion-Assisted Iterative Fluorescence *In Situ* Hybridization). It includes automated image stitching, distributed multi-round image registration, cell segmentation, and distributed spot detection. We also envision this pipeline being adapted for analysis of other image-based spatial transcriptomic data. 
![](/resources/Pipeline.gif)
The pipeline takes advantage of the [n5](https://github.com/saalfeldlab/n5) filesystem to allow for rapid data reading and writing. Fiji-based [n5-viewer](https://github.com/saalfeldlab/n5-viewer) allows for visualization of large image data visualization on a laptop.  

## Pipeline #
An end-to-end analysis [pipeline](https://github.com/JaneliaSciComp/multifish) that takes `czi` image files acquired from Zeiss Z.1 lightsheet microscope and output the transcript spot counts that can be readily used for cell type identification.  

## Modules #

### Stitching #
The previously developed Apache Spark-based high-performance computing pipeline (Gao et al., 2019) was implemented for image stitching. For details on running on local machine, computing cluster or public platforms, please visit  [stitching-spark](https://github.com/saalfeldlab/stitching-spark). 


### Registration #
We developed [BigStream](https://github.com/GFleishman/bigstream) for distributed alignment multi-round FISH data. BigStream first performs fast global affine transformation using a feature-based random sample consensus (RANSAC) algorithm (Fischler and Bolles, 1981). The image volume is then divided into overlapping blocks and another round of feature-based affine transformation was performed, followed by a fast 3D deformable registration [greedypy](https://github.com/GFleishman/greedypy) (Yushkevich, 2016) on each block. 

`bigstream` can be installed with `pip`:\
    `pip install bigstream`

### Segmentation #
[Starfinity](https://github.com/mpicbg-csbd/stardist/tree/refinement) is a deep learning-based automatic 3D segmentation software. Starfinity is an extension of [Stardist](https://github.com/mpicbg-csbd/stardist), an earlier cell detection approach (Schmidt et al., 2018; Weigert et al., 2020) and is based on the dense prediction of cell border distances and their subsequent aggregation into pixel affinities. 

`Starfinity`  can be installed with  `pip`:\
    `pip install git+https://github.com/mpicbg-csbd/stardist@refinement` 

### Spot detection #
We developed hAirlocalize, a modification based on the matlab spot detection algorithm, Airlocalize ([Lionnet et al., 2011](https://www.nature.com/articles/nmeth.1551)) to allow for distributed spot detection in large datasets.  

To use `hAirlocalize`:

## Additional information #

### Post processing
Code used for [Post processing](https://github.com/multiFISH/EASI-FISH/tree/master/post_processing), such as assign spots, filter ROIs, dense spots analysis, intensity measurements, lipofuscin subtraction are included in the repository. 

### Example data #
Example [datasets]() are provided for testing.  

