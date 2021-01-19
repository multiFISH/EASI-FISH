# EASI-FISH analysis toolbox
Table of Contents
=================
   * [Description](#description)
      * [Stitching](#stitching)
      * [Registration](#registration)
      * [Segmentation](#segmentation)
      * [Spot detection](#spot-detection)
      * [Additional analysis](#additional-analysis)     
   * [Installation and Examples](#installation_and_examples)
   * [Pipeline ](#pipeline)
# Description #
<img src="https://github.com/multiFISH/EASI-FISH/blob/master/docs/png/EASI-FISH_example.png" align="right" width="600">
Characterizing the spatial organization and morphological properties of molecular cell types is fundamental for underpinning tissue function. Expansion-Assisted Iterative Fluorescence In Situ Hybridization (EASI-FISH) performs large-scale, multi-round, high resolution FISH imaging that allows gene expression profiling of molecular cell types in thick tissue specimens (300µm). Example dataset from EASI-FISH on right (scale bar: 100µm)




To address the challenge in analyzing multi-terabyte imaging data EASI-FISH produces, we provide a computational pipeline that allows for rapidly processing of such datasets.  The pipeline, which includes automated image stitching, multi-round image registration, cell segmentation, and spot extraction could facilitate adoption of high-plex FISH as a routine laboratory method for tissue analysis. We also envision this pipeline being adapted for analysis of other image-based spatial transcriptomic data. 
![](/docs/png/Pipeline.png)
The pipeline takes advantage of the [n5](https://github.com/saalfeldlab/n5) filesystem to allow for rapid data reading and writing. Fiji-based [n5-viewer](https://github.com/saalfeldlab/n5-viewer) allows for visualization of large image data visualization on a laptop.  
# Stitching #
For large sample volumes, multiple sub-volumes (tiles) need to be sequentially acquired and computational stitched into a single large image. The previously developed Apache Spark-based high-performance computing pipeline [stitching-spark](https://github.com/saalfeldlab/stitching-spark) (Gao et al., 2019) is used for image stitching. The pipeline first performed a flat-field correction for each tile to account for intensity variations and then stitched the intensity-corrected tiles together using an automated and iteratively refined prediction model based on tile coordinates. 

<img src="https://github.com/multiFISH/EASI-FISH/blob/master/docs/png/Stitching.png" width="500">
# Registration #
To register image volumes across multiple rounds of FISH, a robust and fully automatic non-rigid registration pipeline[BigSTREAM](https://github.com/GFleishman/stream) is developed. The analysis pipeline first performs fast global affine transformation using a feature-based random sample consensus (RANSAC) algorithm (Fischler and Bolles, 1981). The image volume is then divided into overlapping blocks and another round of feature-based affine transformation was performed, followed by a fast 3D deformable registration [greedypy](https://github.com/GFleishman/greedypy) (python implementation)(Yushkevich, 2016) on each block. 
<img src="https://github.com/multiFISH/EASI-FISH/blob/master/docs/png/Registration.png" align="left"> 
# Segmentation #
Accurate segmentation of in situ-stained volumetric (3D) fluorescence image data has been a long-standing challenge that can considerably degrade the accuracy of multiplexed FISH analysis pipelines.To overcome this challenge, a deep learning-based automatic 3D segmentation, called [Starfinity](https://github.com/mpicbg-csbd/stardist/tree/refinement) was developed. Starfinity is an extension of [Stardist](https://github.com/mpicbg-csbd/stardist), an earlier cell detection approach (Schmidt et al., 2018; Weigert et al., 2020) and is based on the dense prediction of cell border distances and their subsequent aggregation into pixel affinities. 
# Spot detection #
A spot detection pipeline based on Airlocalize (Lionnet et al., 2011) was developed that allows for parallel processing of chunked overlapping image blocks simultaneously. 
![](/docs/png/Spot_detection.png)

# Additional analysis #
[Additional analysis](https://github.com/multiFISH/EASI-FISH/tree/master/docs/post_processing), such as dense spots analysis, intensity measurements are provided. 

# Installation and Examples #
A example [dataset]() is provided for testing the pipeline. 

# Pipeline # 
We also provide an end-to-end analysis [pipeline](https://github.com/JaneliaSciComp/multifish) in here. 

