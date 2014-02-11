#!/bin/bash

## Brings T1-associated files to single participant / single session space

cd /tmp
for SUB in ${SUBJECTS}; do
    DIR_SESS=`ls -d -- ${DIR_DATA}/${DIR_EXPT}/${SUB}/${DATA_TYPE}/*/`
    for SESS in `basename ${DIR_SESS}`; do
        
        DIR=`echo ${DIR_DATA}/${DIR_EXPT}/${SUB}/${DATA_TYPE}/${SESS}`
        DIR_T1=`echo ${DIR_DATA}/${DIR_EXPT}/${SUB}/T1/${SESS}`
        # register aparc atlas to EPI
        if [ ! -f ${DIR}/anat_aparc_reg.nii.gz ]; then
            flirt -in ${DIR_T1}/anat_aparc_brain.nii.gz \
                  -ref ${DIR}/anat_EPI_brain.nii.gz \
                  -applyxfm -init ${DIR}/mat_T1_to_EPI.mat \
                  -interp nearestneighbour \
                  -out ${DIR}/anat_aparc_reg.nii.gz
        fi

        # register aparc2009 atlas to EPI
        if [ ! -f ${DIR}/anat_aparc2009_reg.nii.gz ]; then
            flirt -in ${DIR_T1}/anat_aparc2009_brain.nii.gz \
                  -ref ${DIR}/anat_EPI_brain.nii.gz \
                  -applyxfm -init ${DIR}/mat_T1_to_EPI.mat \
                  -interp nearestneighbour \
                  -out ${DIR}/anat_aparc2009_reg.nii.gz
        fi
    done
done
cd ${DIR_PIPE}

## JDV Jan 30 2014