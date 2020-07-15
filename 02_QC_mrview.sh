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
echo "Now viewing the FODs (Blue=WM; Green=GM; Red=CSF)"
mrview vf.mif -odf.load_sh wmfod.mif


### Quality checks for Step 3 ###

# Check alignment of the 5 tissue types before and after alignment (new alignment in red, old alignment in blue)
echo "Checking alignment between grey matter alignment before (blue) and after (red)"
mrview dwi_den_preproc_unbiased.mif -overlay.load 5tt_nocoreg.mif -overlay.colourmap 2 -overlay.load 5tt_coreg.mif -overlay.colourmap 1

# Check the seed region (should match up along the GM/WM boundary)
echo "Checking alignment of the seed region with the GM/WM boundary"
mrview dwi_den_preproc_unbiased.mif -overlay.load gmwmSeed_coreg.mif


### Quality checks for Step 4 ###

# View the tracks in mrview
echo "Now viewing the tracks in mrview (red=left-to-right; blue=bottom-to-top; green=forward-to-back)"
mrview dwi_den_preproc_unbiased.mif -tractography.load smallerTracks_200k.tck

# View the sifted tracks in mrview
# Uncomment the following line of code if you used tcksift; otherwise, tcksift2 will output a text file with weightings that are used for later commands (e.g., creating the connectome)
#mrview dwi_den_preproc_unbiased.mif -tractography.load sift_1mio.tck

cd dwifslpreproc-tmp-*
totalSlices=`mrinfo dwi.mif | grep Dimensions | awk '{print $6 * $8}'`
totalOutliers=`awk '{ for(i=1;i<=NF;i++)sum+=$i } END { print sum }' dwi_post_eddy.eddy_outlier_map`
echo "If the following number is greater than 10, you may have to discard this subject because of too much motion or corrupted slices"
echo "scale=5; ($totalOutliers / $totalSlices * 100)/1" | bc | tee percentageOutliers.txt
cd ..
