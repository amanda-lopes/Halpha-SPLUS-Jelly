# Halpha-SPLUS-Jelly

Repository contains the four stages of analysis for Jellyfish galaxies using S-PLUS images.

For such, it uses MontagePy to combine S-PLUS images (r, J0660, i) to create a detection image, SEP to handle the detection objects, masking and estimating the flux radius of the objects, together with splusdata to download the images and splus_datacube to calibrate and create the datacubes.
The analysis is divided in 4 steps:

1) Creates the datacube, detects the sources in the r-image, and masks all sources but the main galaxy.

2) Estimates the radius at which the flux of the galaxies is 20%, 50%, 70% and 90%.

3) Creates the Halpha+[NII] emission line maps applying a Three Filter Method to each pixel of the image.

4) Corrects the integrated emission line flux by extinction and [NII], and converts the Halpha fluxes to star formation rates (SFR).
