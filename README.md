# Halpha-SPLUS-Jelly

Contains the main codes to analyze the Jellyfish galaxies using S-PLUS images.

For such, it uses MontagePy to combine S-PLUS images (r, J0660, i) to create a detection image, SEP to handle the detection objects, masking and estimating the flux radius of the objects, together with splusdata package to download the images and splus_datacube to calibrate and create the datacubes.
The main code to SFR is halpha_sfr_jelly.py.

The analysis is divided in 4 steps:

1) Pre-processing: creates the calibrates the images, creates the datacube, detects the sources in the r-image, and creates masks all sources but the main galaxy. Moreover, performs a PSF-correction by convolving all images with the PSF kernel from the worst PSF image, and applies a low-frequency filter named Butterworth;

2) Estimates the radius at which the flux of the galaxies is 20%, 50%, 70% and 90%. Examples of these radius in respect to the r-band image is located at the folder "jelly_radius_figs";

3) Creates the Halpha+[NII] emission line maps applying a Three Filter Method (Vilella-Rojo et al. 2015) to each pixel of the image. Examples of Halpha+[NII] maps in folder "Halpha_nii_maps";

4) The fluxes in the Halpha+[NII] maps are integrated within 90% radius, estimated at step 2. A second run is made, adding the G-band, in order to correct the Halpha+[NII] emission fluxes by extinction and [NII] applying Eq. (20) and (21) from Vilella-Rojo et al. (2015), respectively. Then, the Halpha fluxes are converted to star formation rates (SFR) using Kennicutt et al. (1998).


Observations: The Halpha fluxes may be over/underestimate mainly due to the dust correction not being the most appropriate. More studies are being done in order to obtain more reliable corrections, in particular related to cases where the HII are not spatially resolved.
