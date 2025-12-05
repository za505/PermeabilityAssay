# PermeabilityAssay
This repository comprises code used to analyze data from the permeability assay in the following papers: "Peptidoglycan turnover promotes active transport of protein through the bacterial cell wall" and "Wall teichoic acids regulate peptidoglycan synthesis by paving cell wall nanostructure"

This repository contains scripts written in MATLAB (2025a). The general workflow uses FIJI and MATLAB to perform image preprocessing, cell tracking, fluorescence measurement, photobleach correction, and fluorescence normalization for the .nd2 movie from each experiment. 
Functions needed to execute each MATLAB script are either defined at the bottom of the .m file or saved as a script file in the 'functions' directory. 

## Preprocessing
First, in FIJI, split the phase and GFP channels. Save the .tiff image stacks in separate directories. Also, crop a portion of the phase movie (which contains features from which it is possible to track drift) and save the reference image stack in its own directory.

Run the 'ImageAlign.m' script. A new directory of aligned imaged stacks will be generated. 
Duplicate the aligned phase directory and rename the ending (i.e. 'erased').

Run the 'EraseImagePart.m' script to delete untrackable cells and remove background (i.e. pillars from the microfluidics device) from the images in the 'erased' directory. 

## Analysis
Track cells from the phase image stack (with background removed) by running BacTrack.m, which generates a '_BT.mat' file as its output. The BT.mat file includes the pixel coorinates for each cell.

Measure the fluorescence intensity from the aligned GFP image stack by running DiffusionMeasure.m, which generate a '_dm.mat' file as its output. The dm.mat file includes the mean intensity value for each cell based on the pixel coordinates in BT.mat. 

Based on your fluorophore--mNeonGreen or FlAsH--use the appropriate DataNormalization.mat script to correct for photobleaching and normalize fluorescence. 

After fluorescence normalization, the '_nct.mat' output file from a set of experiments can be compiled and stored as a separate .mat file for easier data analysis.

FinalPlots.m is the code used the generate the plots shown in Fig 3A-C of the paper, "Wall teichoic acids regulate peptidoglycan synthesis by paving cell wall nanostructure"


