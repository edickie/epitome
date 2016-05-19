## fs2hcp of a glmdir
## this build from some of fs2hcp to try to convert freesurfer glm outputs for visulization with hcp

## template_hcp_sub

hcp_folder=/scratch/edickie/fs_glm_to_hcp/hcp
fs_folder=/scratch/edickie/fs_glm_to_hcp/freesurfer
Subject=cvs_avg35_inMNI152
fs_glm_dir=/projects/laura/current/Freesurfer/GLM_Age+Gender/rh.glm_CT_TASIT3_AG.glmdir/group-x-gender
ana_name="CT_TASIT3_AG"

# info about output files taken from https://surfer.nmr.mgh.harvard.edu/fswiki/FsTutorial/GroupAnalysis
# uncorrected values
# C.dat -- original contrast matrix (text)
# cnr.mgh -- contrast-to-noise ratio
# efficiency.mgh -- statistical efficiency for the contrast
# F.mgh -- F ratio of contrast  (surface overlay)
# gamma.mgh -- contrast effect size (surface overlay)
# gammvar.mgh --contrast variance
# maxvox.dat -- voxel with the maximum statistic
# pcc.mgh -- partial (pearson) correlation coefficien
# sig.mgh -- significance, -log10(pvalue), uncorrected (surface overlay)
# z.mgh -- z-stat that corresponds to the significance
#
# cluster (montecarlo) corrected maps
# cache.th40.neg.sig.pdf.dat -- probability distribution function of clusterwise correction
# cache.th40.neg.sig.cluster.mgh -- cluster-wise corrected map (overlay)
# cache.th40.neg.sig.cluster.summary -- summary of clusters (text)
# cache.th40.neg.sig.masked.mgh -- uncorrected sig values masked by the clusters that survive correction
# cache.th40.neg.sig.ocn.annot -- output cluster number (annotation of clusters)
# cache.th40.neg.sig.ocn.mgh -- output cluster number (segmentation showing where each numbered cluster is)
# cache.th40.neg.sig.voxel.mgh -- voxel-wise map corrected for multiple comparisons at a voxel (rather than cluster) level
# cache.th40.neg.sig.y.ocn.dat -- the average value of each subject in each cluster

sig.mgh
cache.th13.abs.sig.cluster.mgh
cache.th13.abs.sig.ocn.annot
cache.th13.abs.sig.masked.mgh

hemisphere='r'
Hemisphere='R'
FreeSurferFolder=${fs_folder}/fsaverage
AtlasSpaceFolder=${hcp_folder}/${Subject}/MNINonLinear
NativeFolder="Native"


mris_convert -c "$fs_glm_dir"/sig.mgh \
  "$FreeSurferFolder"/surf/"$hemisphere"h.white \
  "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$ana_name".sig.native.shape.gii
wb_command -set-structure "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$wbname".native.shape.gii ${Structure}
wb_command -metric-math "var * -1" "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$wbname".native.shape.gii -var var "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$wbname".native.shape.gii
wb_command -set-map-names "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$wbname".native.shape.gii -map 1 "$Subject"_"$Hemisphere"_"$mapname"
wb_command -metric-palette "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere"."$wbname".native.shape.gii MODE_AUTO_SCALE_PERCENTAGE -pos-percent 2 98 -palette-name Gray_Interp -disp-pos true -disp-neg true -disp-zero true
