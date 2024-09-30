Comprehensive README and code commenting to come. 


# Run CRISP:
This folder contains the files to actually run CRISP on a users images.

## Mask_refinement
This function takes a time course of pixels with masks differentiating
cells and removes pixels within that mask that do not correlate with the rest of the pixels (e.g. likely not a part of the cell.)

Inputs: 
- Images is the image file in matrix form (x,y,t). Currently only accepting 3D matrices (e.g. no z stacks) so z stacks should be fed in separately. 

- Cell Mask is a pixel x pixel matrix of the imaging file. Any pixel that is considered to be a part of a cell has the index of that corresponding cell. For example, if we are looking at a 5x5 pixel images with 1 cell inside, cell mask would be:  [0, 0, 0, 0, 0; 0, 1, 0, 0, 0; 0, 1, 1, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0]

- Opts - Options structure containing thresholds and figure options: 
-- Opt.fig = 0 or 1 for no figures or figures respectively
-- Opt.st_thr = double that defines how many standard deviations to count as bad pixels. If Opt.Thr == 'corr' then Opt.st_thr is the actual threshold for the correlation coefficient. 
-- Opt.Thr = 'corr', 'st', if you want to use a fixed correlation threshold (corr) or a threshold based on the cell's correlation distribution (st). (Manuscript used st to be more generalizable)

## CRISP_annulus_corr
This function takes a nucleus location and itteratively pulls more pixels around that location until the pixels on the outermost edge appear to not be apart of the cell. Whether or not the pixels are a part of the cell is based on the assumption that pixels from the same cell will have extremely similar behavior (minus noise). Therefore, once the behavior of pixels on the outermost radius of the cell deviates too much (based on the threshold (thr), the cell radius and corresponding timecourse is exported.

Inputs:
- images is the image file in matrix form (x,y,t). Currently only accepting. 3D matrices (e.g. no z stacks) so z stacks should be fed in separately.

- NucLoc is the X,Y position of the center of the nucleus (note, if this changes overtime, take the median for the most accurate result).

- opts - Options structure containing thresholds and figure options: 
-- opt.figs = 0 or 1 for no figures or figures respectively
-- opts.score_thr = score threshold
-- opts.th_pix = pixel correlation threshold



# Training and Testing: 

This folder contains code for training and testing CRISP, including masking calcium images, creating ROC curves, and running network analysis. 

In *MaskRefinement_Training.m*, *NetworkAnalysis* function is called -- this function can be downloaded at: https://github.com/jenniferkbriggs/Islet_Heterogeneity

# Plotting:
Various code for plotting results. 


