#!/bin/bash

# Written by Andrew Jahn, University of Michigan, 02.25.2019
# Based on Marlene Tahedl's BATMAN tutorial (http://www.miccai.org/edu/finalists/BATMAN_trimmed_tutorial.pdf)

display_usage() {
	echo "$(basename $0) [Raw Diffusion] [Fieldmap] [AP bvec] [AP bval] [PA bvec] [PA bval] [Anatomical]"
	echo "This script uses MRtrix to analyze diffusion data. It requires 7 arguments: 
		1) The raw diffusion image;
		2) The fieldmap image;
		3) The bvec file for the data acquired in the AP direction;
		4) The bval file for the data acquired in the AP direction;
		5) The bvec file for the data acquired in the PA direction;
		6) The bval file for the data acquired in the PA direction;
		7) The anatomical image"
	}

	if [ $# -le 6 ]
	then
		display_usage
		exit 1
	fi

RAW_DWI=$1
FIELDMAP=$2
BVEC=$3
BVAL=$4
PA_BVEC=$5
PA_BVAL=$6
ANAT=$7


########################### STEP 1 ###################################
# Convert data to .mif format, denoise, and perform Gibbs de-ringing #
######################################################################

mrconvert $RAW_DWI raw_dwi.mif -fslgrad $BVEC $BVAL
dwidenoise raw_dwi.mif dwi_den.mif -noise noise.mif
mrcalc raw_dwi.mif dwi_den.mif -subtract residual.mif
mrdegibbs dwi_den.mif dwi_den_unr.mif -axes 0,1

# Extract the b0 images from the diffusion data acquired in the AP direction
dwiextract dwi_den_unr.mif - -bzero | mrmath - mean mean_b0_AP.mif -axis 3

# Extracts the b0 images for diffusion data acquired in the PA direction
# In the fieldmap image, there are 8 images in this order: AP, PA, AP, PA, AP, PA, AP, PA. Indexing starts at 0
mrconvert $FIELDMAP -coord 3 1 PA.mif
mrconvert PA.mif -fslgrad $PA_BVEC $PA_BVAL - | mrmath - mean mean_b0_PA.mif -axis 3

# Concatenates the b0 images from AP and PA directions to create a paired b0 image
mrcat mean_b0_AP.mif mean_b0_PA.mif -axis 3 b0_pair.mif

# Runs the dwipreproc command, which is a wrapper for eddy and topup. This step takes about 2 hours
dwipreproc dwi_den_unr.mif dwi_den_unr_preproc.mif -pe_dir AP -rpe_pair -se_epi b0_pair.mif -eddy_options " --slm=linear --data_is_shelled"

# Performs bias field correction. Needs ANTs to be installed in order to use the "ants" option (use "fsl" otherwise)
dwibiascorrect -ants dwi_den_unr_preproc.mif dwi_den_unr_preproc_unbiased.mif -bias bias.mif

# Create a mask for future processing steps
dwi2mask dwi_den_unr_preproc_unbiased.mif mask.mif


########################### STEP 2 ###################################
#             Basis function for each tissue type                    #
######################################################################

# Create a basis function from the subject's DWI data. The "dhollander" function is best used for multi-shell acquisitions; it will estimate different basis functions for each tissue type. For single-shell acquisition, use the "tournier" function instead
dwi2response dhollander dwi_den_unr_preproc_unbiased.mif wm.txt gm.txt csf.txt -voxels voxels.mif

# Performs multishell-multitissue constrained spherical deconvolution, using the basis functions estimated above
dwi2fod msmt_csd dwi_den_unr_preproc_unbiased.mif -mask mask.mif wm.txt wmfod.mif gm.txt gmfod.mif csf.txt csffod.mif

# Creates an image of the fiber orientation densities overlaid onto the estimated tissues (Blue=WM; Green=GM; Red=CSF)
# You should see FOD's mostly within the white matter. These can be viewed later with the command "mrview vf.mif -odf.load_sh wmfod.mif"
mrconvert -coord 3 0 wmfod.mif - | mrcat csffod.mif gmfod.mif - vf.mif

# Now normalize the FODs to enable comparison between subjects
mtnormalise wmfod.mif wmfod_norm.mif gmfod.mif gmfod_norm.mif csffod.mif csffod_norm.mif -mask mask.mif


########################### STEP 3 ###################################
#            Create a GM/WM boundary for seed analysis               #
######################################################################

# Convert the anatomical image to .mif format, and then extract all five tissue catagories (1=GM; 2=Subcortical GM; 3=WM; 4=CSF; 5=Pathological tissue)
mrconvert $ANAT anat.mif
5ttgen fsl anat.mif 5tt_nocoreg.mif

# The following series of commands will take the average of the b0 images (which have the best contrast), convert them and the 5tt image to NIFTI format, and use it for coregistration.
dwiextract dwi_den_unr_preproc_unbiased.mif - -bzero | mrmath - mean mean_b0_processed.mif -axis 3
mrconvert mean_b0_processed.mif mean_b0_processed.nii.gz
mrconvert 5tt_nocoreg.mif 5tt_nocoreg.nii.gz

# Uses FSL commands fslroi and flirt to create a transformation matrix for regisitration between the tissue map and the b0 images
fslroi 5tt_nocoreg.nii.gz 5tt_vol0.nii.gz 0 1 #Extract the first volume of the 5tt dataset (since flirt can only use 3D images, not 4D images)
flirt -in mean_b0_processed.nii.gz -ref 5tt_vol0.nii.gz -interp nearestneighbour -dof 6 -omat diff2struct_fsl.mat
transformconvert diff2struct_fsl.mat mean_b0_processed.nii.gz 5tt_nocoreg.nii.gz flirt_import diff2struct_mrtrix.txt
mrtransform 5tt_nocoreg.mif -linear diff2struct_mrtrix.txt -inverse 5tt_coreg.mif

#Create a seed region along the GM/WM boundary
5tt2gmwmi 5tt_coreg.mif gmwmSeed_coreg.mif

########################### STEP 4 ###################################
#                 Run the streamline analysis                        #
######################################################################

# Create streamlines
# Note that the "right" number of streamlines is still up for debate. Last I read from the MRtrix documentation,
# They recommend about 100 million tracks. Here I use 10 million, if only to save time. Read their papers and then make a decision
tckgen -act 5tt_coreg.mif -backtrack -seed_gmwmi gmwmSeed_coreg.mif -select 10000k wmfod_norm.mif tracks_10mio.tck

# Extract a subset of tracks (here, 200 thousand) for ease of visualization
tckedit tracks_10mio.tck -number 200k smallerTracks_200k.tck

# Reduce the number of streamlines with tcksift
tcksift -act 5tt_coreg.mif -term_number 1000000 tracks_10mio.tck wmfod_norm.mif sift_1mio.tck

