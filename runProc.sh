#!/bin/bash

#Based on the fixel-based analysis steps outlined on the MRtrix website (www.mrtrix.org).
#Note that this is a rough draft, and is designed to be used with the BTC_Preop data on openneuro.org. This script analyzes six subjects, three from each group.

for i in CON01 CON02 CON03 PAT01 PAT02 PAT03; do
mkdir $i;
mv *$i*.nii.gz *$i*.bval *$i*.bvec $i;
for_each * : mrconvert IN/*$i*.nii.gz IN/$i_dwi.mif -fslgrad IN/*$i*.bvec IN/*$i*.bval
done

for_each * : mrconvert -force IN/*.nii.gz IN/NAME.mif -fslgrad IN/*.bvec IN/*.bval
for_each * : dwidenoise IN/*.mif IN/dwi_denoised.mif
for_each * : dwifslpreproc -nthreads 4 IN/dwi_denoised.mif IN/dwi_denoised_preproc.mif -rpe_none -pe_dir AP
for_each * : dwibiascorrect ants IN/dwi_denoised_preproc.mif IN/dwi_denoised_preproc_unbiased.mif
for_each * : dwi2response dhollander IN/dwi_denoised_preproc_unbiased.mif IN/response_wm.txt IN/response_gm.txt IN/response_csf.txt
responsemean */response_wm.txt ../group_average_response_wm.txt
responsemean */response_gm.txt ../group_average_response_gm.txt
responsemean */response_csf.txt ../group_average_response_csf.txt
for_each * : mrgrid IN/dwi_denoised_preproc_unbiased.mif regrid -vox 1.25 IN/dwi_denoised_preproc_unbiased_upsampled.mif
for_each * : dwi2mask IN/dwi_denoised_preproc_unbiased_upsampled.mif IN/dwi_mask_upsampled.mif
for_each * : dwi2fod msmt_csd IN/dwi_denoised_preproc_unbiased_upsampled.mif ../group_average_response_wm.txt IN/wmfod.mif ../group_average_response_gm.txt IN/gm.mif  ../group_average_response_csf.txt IN/csf.mif -mask IN/dwi_mask_upsampled.mif
for_each * : mtnormalise IN/wmfod.mif IN/wmfod_norm.mif IN/gm.mif IN/gm_norm.mif IN/csf.mif IN/csf_norm.mif -mask IN/dwi_mask_upsampled.mif
mkdir -p ../template/fod_input
mkdir ../template/mask_input
for_each * : ln IN/wmfod_norm.mif ../template/fod_input/PRE.mif
for_each * : ln IN/dwi_mask_upsampled.mif ../template/mask_input/PRE.mif
#for_each `ls -d P* | sort -R | tail -20` : ln -s IN/wmfod_norm.mif ../template/fod_input/PRE.mif ";" ln -s IN/dwi_mask_upsampled.mif ../template/mask_input/PRE.mif
#for_each `ls -d C* | sort -R | tail -20` : ln -s IN/wmfod_norm.mif ../template/fod_input/PRE.mif ";" ln -s IN/dwi_mask_upsampled.mif ../template/mask_input/PRE.mif
population_template ../template/fod_input -mask_dir ../template/mask_input ../template/wmfod_template.mif -voxel_size 1.25
for_each * : mrregister IN/wmfod_norm.mif -mask1 IN/dwi_mask_upsampled.mif ../template/wmfod_template.mif -nl_warp IN/subject2template_warp.mif IN/template2subject_warp.mif
for_each * : mrtransform IN/dwi_mask_upsampled.mif -warp IN/subject2template_warp.mif -interp nearest -datatype bit IN/dwi_mask_in_template_space.mif
mrmath */dwi_mask_in_template_space.mif min ../template/template_mask.mif -datatype bit
fod2fixel -mask ../template/template_mask.mif -fmls_peak_value 0.06 ../template/wmfod_template.mif ../template/fixel_mask
for_each * : mrtransform IN/wmfod_norm.mif -warp IN/subject2template_warp.mif -reorient_fod no IN/fod_in_template_space_NOT_REORIENTED.mif
for_each * : fod2fixel -mask ../template/template_mask.mif IN/fod_in_template_space_NOT_REORIENTED.mif IN/fixel_in_template_space_NOT_REORIENTED -afd ed.mif
for_each * : fixelreorient IN/fixel_in_template_space_NOT_REORIENTED IN/subject2template_warp.mif IN/fixel_in_template_space
for_each * : fixelcorrespondence IN/fixel_in_template_space/fd.mif ../template/fixel_mask ../template/fd PRE.mif
for_each * : warp2metric IN/subject2template_warp.mif -fc ../template/fixel_mask ../template/fc IN.mif
mkdir ../template/log_fc
cp ../template/fc/index.mif ../template/fc/directions.mif ../template/log_fc
for_each * : mrcalc ../template/fc/IN.mif -log ../template/log_fc/IN.mif
mkdir ../template/fdc
cp ../template/fc/index.mif ../template/fdc
cp ../template/fc/directions.mif ../template/fdc
for_each * : mrcalc ../template/fd/IN.mif ../template/fc/IN.mif -mult ../template/fdc/IN.mif
cd ../template
tckgen -angle 22.5 -maxlen 250 -minlen 10 -power 1.0 wmfod_template.mif -seed_image template_mask.mif -mask template_mask.mif -select 2000000 -cutoff 0.06 tracks_20_million.tck
tcksift tracks_20_million.tck wmfod_template.mif tracks_2_million_sift.tck -term_number 200000
fixelconnectivity fixel_mask/ tracks_2_million_sift.tck matrix/
fixelfilter fd smooth fd_smooth -matrix matrix/
fixelfilter log_fc smooth log_fc_smooth -matrix matrix/
fixelfilter fdc smooth fdc_smooth -matrix matrix/
fixelcfestats fd_smooth/ files.txt design_matrix.txt contrast_matrix.txt matrix/ stats_fd/
