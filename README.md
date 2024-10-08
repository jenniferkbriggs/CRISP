#CRISP: CRISP: Correlation-Refined Image Segmentation Process

This repository contains all code used in the development and testing of CRISP. The code is sorted into folders: 
1. [Run CRISP](#Run_CRISP)
2. [Training and Testing](#Training_and_Testing)
3. [Plotting](#Plotting)

*Manuscript Abstract*
Calcium imaging enables real-time recording of cellular activity across various biological contexts. To assess the activity of individual cells, scientists typically manually outline the cells based on visual inspection. This manual cell masking introduces potential user error. To ameliorate this error, we developed the Correlation-Refined Image Segmentation Process (CRISP), a two-part automated algorithm designed to both enhance the accuracy of user-drawn cell masks and to automatically identify cell masks. We developed and tested CRISP on calcium images of densely packed β-cells within the islet of Langerhans. Because these β-cells are densely packed within the islet, traditional clustering-based image segmentation methods struggle to identify individual cell outlines. Using β-cells isolated from two different mouse phenotypes and imaged on two different confocal microscopes, we show that CRISP is generalizable and accurate. To test the benefit of using CRISP in functional biological analyses, we show that CRISP improves accuracy of functional network analysis and utilize CRISP to characterize the distribution of β-cell size.


## Run_CRISP:
This folder contains the files to run CRISP on a user's images.

### Mask_refinement
This function takes a time course of pixels with masks differentiating
cells and removes pixels within that mask that do not correlate with the rest of the pixels (e.g. likely not a part of the cell.)

Inputs: 
- Images is the image file in matrix form (x,y,t). Currently only accepting 3D matrices (e.g. no z stacks) so z stacks should be fed in separately. 

- Cell Mask is a pixel x pixel matrix of the imaging file. Any pixel that is considered to be a part of a cell has the index of that corresponding cell. For example, if we are looking at a 5x5 pixel images with 1 cell inside, cell mask would be:  [0, 0, 0, 0, 0; 0, 1, 0, 0, 0; 0, 1, 1, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0]

- Opts - Options structure containing thresholds and figure options:
    - Opt.fig = 0 or 1 for no figures or figures respectively
    - Opt.st_thr = double that defines how many standard deviations to count as bad pixels. If Opt.Thr == 'corr' then Opt.st_thr is the actual threshold for the correlation coefficient.
    - Opt.Thr = 'corr', 'st', if you want to use a fixed correlation threshold (corr) or a threshold based on the cell's correlation distribution (st). (Manuscript used st to be more generalizable)

### CRISP_annulus_corr
This function takes a nucleus location and iteratively pulls more pixels around that location until the pixels on the outermost edge appear to not be a part of the cell. Whether or not the pixels are a part of the cell is based on the assumption that pixels from the same cell will have extremely similar behavior (minus noise). Therefore, once the behavior of pixels on the outermost radius of the cell deviates too much (based on the threshold (thr), the cell radius and corresponding timecourse are exported.

Inputs:
- images is the image file in matrix form (x,y,t). Currently only accepting. 3D matrices (e.g. no z stacks) so z stacks should be fed in separately.

- NucLoc is the X,Y position of the center of the nucleus (note, if this changes over time, take the median for the most accurate result).

- opts - Options structure containing thresholds and figure options:
    - opt.figs = 0 or 1 for no figures or figures respectively
    - opts.score_thr = score threshold
    - opts.th_pix = pixel correlation threshold



## Training_and_Testing: 

This folder contains code for training and testing CRISP, including masking calcium images, creating ROC curves, and running network analysis. 

In *MaskRefinement_Training.m*, *NetworkAnalysis* function is called -- this function can be downloaded at: https://github.com/jenniferkbriggs/Islet_Heterogeneity

## Plotting:
Various code for plotting results. 


