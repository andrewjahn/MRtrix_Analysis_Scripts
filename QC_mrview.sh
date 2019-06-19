#!/bin/bash

# These commands are for quality-checking your diffusion data


### Quality checks for Step 2 ###

# Views the voxels used for FOD estimation
echo "Now viewing the voxels used for FOD estimation (Blue=WM; Green=GM; Red=CSF)"
mrview dwi_den_preproc_unbiased.mif -overlay.load voxels.mif

# Views the response functions for each tissue type. The WM function should flatten out at higher b-values, while the other tissues should remain spherical
echo "Now viewing response function for white matter (press right arrow key to view response function for different shells)"
shview wm.txt
echo "Now viewing response function for grey matter"
shview gm.txt
echo "Now viewing response function for CSF"
shview csf.txt

# Views the FODs overlaid on the tissue types (Blue=WM; Green=GM; Red=CSF)
mrview vf.mif -odf.load_sh wmfod.mif


### Quality checks for Step 3 ###

# Check alignment of the 5 tissue types before and after alignment (new alignment in red, old alignment in blue)
mrview dwi_den_preproc_unbiased.mif -overlay.load 5tt_nocoreg.mif -overlay.colourmap 2 -overlay.load 5tt_coreg.mif -overlay.colourmap 1

# Check the seed region (should match up along the GM/WM boundary)
mrview dwi_den_preproc_unbiased.mif -overlay.load gmwmSeed_coreg.mif


### Quality checks for Step 4 ###

# View the tracks in mrview
mrview dwi_den_preproc_unbiased.mif -tractography.load smallerTracks_200k.tck

# View the sifted tracks in mrview
# Uncomment the following line of code if you used tcksift; otherwise, tcksift2 will output a text file with weightings that are used for later commands (e.g., creating the connectome)
#mrview dwi_den_preproc_unbiased.mif -tractography.load sift_1mio.tck

