#!/bin/bash


#Convert the labels of the FreeSurfer parcellation to a format that MRtrix understands. This requires recon-all to have been run on the subject
labelconvert -force sub-01_MRtrix/mri/aparc+aseg.mgz $FREESURFER_HOME/FreeSurferColorLUT.txt ~/mrtrix3/share/mrtrix3/labelconvert/fs_default.txt sub-01_parcels.mif

# Add a line here for mrtransform to register the parcels to the diffusion data?

# mrtransform sub-01_parcels.mif -interp nearest –linear diff2struct_mrtrix.txt –inverse –datatype uint32 sub-01_parcels_coreg.mif

#Create a whole-brain connectome, representing the streamlines between each parcellation pair in the atlas (in this case, 84x84). The "symmetric" option will make the lower diagonal the same as the upper diagonal, and the "scale_invnodevol" option will scale the connectome by the inverse of the size of the node 
#tck2connectome -symmetric -zero_diagonal -scale_invnodevol -tck_weights_in sift_1M.txt sub-01_parcels.mif sub-01_parcels.csv -out_assignment assignments_sub-01_parcels.csv
tck2connectome -symmetric -zero_diagonal -scale_invnodevol -tck_weights_in sift_1M.txt tracks_10M.tck sub-01_parcels.mif sub-01_parcels.csv -out_assignment assignments_sub-01_parcels.csv

#Creates a tract file between the specified nodes that can then be visualized in mrview. Replace the "8,10" pair after the "nodes" option with the labels in ~/mrtrix3/share/mrtrix3/labelconvert/fs_default.txt that you are interested in
connectome2tck -nodes 8,10 -exclusive sift_1mio.tck assignments_sub-01_parcels.csv test
