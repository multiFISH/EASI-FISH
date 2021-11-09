# EASI-FISH analysis toolbox # 
[![DOI](https://zenodo.org/badge/319982517.svg)](https://zenodo.org/badge/latestdoi/319982517)

Expansion-Assisted Iterative-FISH defines lateral hypothalamus spatio-molecular organization
Yuhan Wang, Mark Eddison, Greg Fleishman, Martin Weigert, Shengjin Xu, Fredrick E. Henry, Tim Wang, Andrew L. Lemire, Uwe Schmidt, Hui Yang, Konrad Rokicki, Cristian Goina, Karel Svoboda, Eugene W. Myers, Stephan Saalfeld, Wyatt Korff, Scott M. Sternson, Paul W. Tillberg &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
bioRxiv 2021.03.08.434304; doi: https://doi.org/10.1101/2021.03.08.434304

See [here](https://github.com/JaneliaSciComp/multifish) for active EASI-FISH pipeline updates. 
## Table of Contents #
   * [Description](#description)
   * [Pipeline](#pipeline)
   * [Modules](#modules)
      * [Stitching](#stitching)
      * [Registration](#registration)
      * [Segmentation](#segmentation)
      * [Spot detection](#spot-detection)
   * [Additional information](#additional-information)
      * [Visualization](#visualization)
      * [Post processing](#post-processing)  
      * [Example data](#example-data)

## Description #
This workflow is used to analyze large-scale, multi-round, high-resolution image data acquired using EASI-FISH (Expansion-Assisted Iterative Fluorescence *In Situ* Hybridization). It takes advantage of the [n5](https://github.com/saalfeldlab/n5) filesystem to allow for rapid and parallel data reading and writing. It performs automated image stitching, distributed and highly accurate multi-round image registration, 3D cell segmentation, and distributed spot detection. We also envision this workflow being adapted for analysis of other image-based spatial transcriptomic data. 
![](/resources/Pipeline.gif)
 
## Pipeline #
We build a self-contained, highly flexible, and platform agnostic computational [pipeline](https://github.com/JaneliaSciComp/multifish), which supports turnkey EASI-FISH data analysis on local machines and the LSF compute cluster. The pipeline is freely available, open source, and modular. It can rapidly process large datasets greater than 10 TB in size with minimal manual intervention. The pipeline can be used to analyze EASI-FISH dataset end-to-end. It takes `czi` image files acquired from Zeiss Z.1 lightsheet microscope as input and outputs 1) processed image data at different scales and 2) transcript counts that can be readily used for cell type identification. The pipeline also provides options to run individual analysis modules, such as image stitching or registration. 

## Modules #

### Stitching #
For imaging large volumes, multiple sub-volumes (tiles) are sequentially acquired, followed by computational stitching into a single large image. We used an Apache Spark-based high-performance stitching pipeline ([Gao et al., 2019](https://science.sciencemag.org/content/363/6424/eaau8302.long)). The pipeline automatically performs a flat-field correction for each tile to account for intensity variations across the lightsheet. It then derives the globally optimal translation for each tile that minimizes the sum of square distances to competing optimal pairwise translations estimated by phase-correlation ([Preibisch et al., 2009](https://academic.oup.com/bioinformatics/article/25/11/1463/332497)). The stitching module can be executed with the EASI-FISH pipeline (see above). For additional details, please see [stitching-spark](https://github.com/saalfeldlab/stitching-spark). 


### Registration #
We developed [BigStream](https://github.com/GFleishman/bigstream) for robust and fully automated non-rigid registration of multi-round FISH data. BigStream first performs fast global affine transformation using a feature-based random sample consensus (RANSAC) algorithm. The image volume is then divided into overlapping blocks and another round of feature-based affine transformation is performed, followed by a fast 3D deformable registration [greedypy](https://github.com/GFleishman/greedypy) ([Yushkevich, 2016](https://github.com/pyushkevich/greedy)) on each block. Bigstream can be executed as part of the EASI-FISH pipeline. It also can be installed and used seperately. 

`bigstream` can be installed with `pip`:
```
   pip install bigstream
```
### Segmentation #
[Starfinity](https://github.com/mpicbg-csbd/stardist/tree/refinement) is a deep learning-based automatic 3D segmentation software. Starfinity is an extension of [Stardist](https://github.com/mpicbg-csbd/stardist), an earlier cell detection approach (Schmidt et al., 2018; Weigert et al., 2020) and is based on the dense prediction of cell border distances and their subsequent aggregation into pixel affinities. A starfinity [model](https://doi.org/10.25378/janelia.13624268) was trained to predict cell body shapes from DAPI-stained RNA images and is provided for testing. Starfinity can be executed as part of the EASI-FISH pipeline. It can also be installed and used independently. 

`Starfinity`  can be installed with  `pip`:
```
   pip install git+https://github.com/mpicbg-csbd/stardist@refinement
```
For training new Starfinity models, the [augmend](https://github.com/stardist/augmend) and [gputools](https://github.com/maweigert/gputools)(optional) packages need to be installed.

### Spot detection #
We developed hAirlocalize, a distributed spot detection method based on the MATLAB spot detection algorithm, [Airlocalize](https://github.com/timotheelionnet/AIRLOCALIZE)([Lionnet et al., 2011](https://www.nature.com/articles/nmeth.1551)) to allow rapid spot detection on full-resolution large image datasets. hAirlocalize can be executed independently or as part of the EASI-FISH pipeline (see above). For independent execution, we recommend working with the n5 filesystem due to large file size.

## Additional information #

### Visualization #
Fiji-based [n5-viewer](https://github.com/saalfeldlab/n5-viewer) can be used for large image dataset visualization on local machines. The workflow also outputs processed intermediate image data in the stitching (`n5`), registration (`n5`) and segmentation (`tif`) steps. For inspection of spot extracted with hAirlocalize, we recommend the python-based multi-dimensional image viewer, [napari](https://napari.org/). Example [notebooks](https://github.com/multiFISH/EASI-FISH/tree/master/data_visualization) are provided. 

### Post processing
Code used for [Post processing](https://github.com/multiFISH/EASI-FISH/tree/master/data_processing), such as assign spots, cell morphological measurements, dense spot analysis, FISH signal intensity measurements, lipofuscin subtraction are included in this repository. 

### Example data #
EASI-FISH [example datasets](https://doi.org/10.25378/janelia.c.5276708.v1) are provided for software testing. For instructions on performing a demo run with the example data using the end-to-end analysis pipeline, see [here](https://github.com/JaneliaSciComp/multifish). 

