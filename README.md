# Halpha-SPLUS-Jelly

Contains the main codes to analyze the Jellyfish galaxies using S-PLUS images.

For such, it uses MontagePy to combine S-PLUS images (r, J0660, i) to create a detection image, SEP to handle the detection objects, masking and estimating the flux radius of the objects, together with splusdata package to download the images and splus_datacube to calibrate and create the datacubes.
The main code to obtain the SFR is halpha_sfr_jelly.py. This code includes all the necessary steps described below.

The analysis is divided in 4 steps:

1) Pre-processing: creates the calibrates the images, creates the datacube, detects the sources in the r-image, and creates masks all sources but the main galaxy. Moreover, performs a PSF-correction by convolving all images with the PSF kernel from the worst PSF image, and applies a low-frequency filter named Butterworth;

2) Estimates the radius at which the flux of the galaxies is 20%, 50%, 70% and 90%. Examples of these radius in respect to the r-band image is located at the folder "jelly_radius_figs";

3) Creates the Halpha+[NII] emission line maps applying a Three Filter Method (Vilella-Rojo et al. 2015) to each pixel of the image. Remove pixels with larger errors (Halpha+[NII] flux/ Halpha+[NII] flux_err) <2. Examples of Halpha+[NII] maps in folder "Halpha_nii_maps";

4) The fluxes in the Halpha+[NII] maps are integrated within 90% radius, estimated at step 2. Then, the emission fluxes are converted to star formation rates (SFR) using Kennicutt et al. (1998). In order to correct by extinction and [NII], the final SFR is estimated following Eq.(19) from Kouroumpatzakis et al. (2021). 


Observations: The SFR may be underestimate mainly due to the dust correction not being the most appropriate. More studies are being done in order to obtain more reliable corrections, in particular related to cases where the HII regions are not spatially resolved.
