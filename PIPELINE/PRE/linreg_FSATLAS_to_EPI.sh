#!/bin/bash

## Brings T1-associated files to single participant / single session space

cd /tmp
for SUB in ${SUBJECTS}; do
    DIR_SESS=`ls -d -- ${DIR_DATA}/${DIR_EXPT}/${SUB}/${DATA_TYPE}/*/`
    for SESS in ${DIR_SESS}; do
        SESS=`basename ${SESS}`        
        DIR=`echo ${DIR_DATA}/${DIR_EXPT}/${SUB}/${DATA_TYPE}/${SESS}`
        DIR_T1=`echo ${DIR_DATA}/${DIR_EXPT}/${SUB}/T1/${SESS}`
        # register aparc atlas to EPI
        if [ ! -f ${DIR}/anat_aparc_reg.nii.gz ]; then
            3dAllineate \
                -prefix ${DIR}/anat_aparc_reg.nii.gz \
                -input ${DIR_T1}/anat_aparc_brain.nii.gz \
                -1Dmatrix_apply ${DIR}/mat_T1_to_EPI.aff12.1D \
                -master ${DIR}/anat_EPI_mask.nii.gz \
                -float \
                -interp NN \
                -final NN
        fi

        # register aparc2009 atlas to EPI
        if [ ! -f ${DIR}/anat_aparc2009_reg.nii.gz ]; then
            3dAllineate \
                -prefix ${DIR}/anat_aparc2009_reg.nii.gz \
                -input ${DIR_T1}/anat_aparc2009_brain.nii.gz \
                -1Dmatrix_apply ${DIR}/mat_T1_to_EPI.aff12.1D \
                -master ${DIR}/anat_EPI_mask.nii.gz \
                -float \
                -interp NN \
                -final NN
        fi
    done
done
cd ${DIR_PIPE}
