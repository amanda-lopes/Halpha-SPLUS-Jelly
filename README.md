# Halpha-SPLUS-Jelly

Repository contains main codes to analyze the Jellyfish galaxies using S-PLUS images.

For such, it uses MontagePy to combine S-PLUS images (r, J0660, i) to create a detection image, SEP to handle the detection objects, masking and estimating the flux radius of the objects, together with splusdata package to download the images and splus_datacube to calibrate and create the datacubes.
The analysis is divided in 4 steps:

1) Pre-processing: creates the calibrates the images, creates the datacube, detects the sources in the r-image, and creates masks all sources but the main galaxy. Moreover, performs a PSF-correction by convolving all images with the PSF kernel from the worst PSF image, and applies a low-frequency filter named Butterworth;

2) Estimates the radius at which the flux of the galaxies is 20%, 50%, 70% and 90%. Examples of these radius in respect to the r-band image is located at the folder "jelly_radius_figs";

3) Creates the Halpha+[NII] emission line maps applying a Three Filter Method (Vilella-Rojo et al. 2015) to each pixel of the image. Examples of Halpha+[NII] maps in folder "Halpha_nii_maps";

4) Corrects the integrated emission line flux by extinction and [NII], and converts the Halpha fluxes to star formation rates (SFR).
