#!/bin/bash

## Create Regressors Module: Creates a series of regressors from fMRI data and 
## a freesurfer segmentation

# white matter + eroded mask
# ventricles + eroded mask
# grey matter mask
# brain stem mask
# dia brain mask
# draining vessels mask
# local white matter regressors + lagged
# ventricle regressors + lagged
# draining vessel regressors + lagged

cd /tmp
for SUB in ${SUBJECTS}; do
    DIR_SESS=`ls -d -- ${DIR_DATA}/${DIR_EXPT}/${SUB}/${DATA_TYPE}/*/`
    for SESS in ${DIR_SESS}; do
        SESS=`basename ${SESS}`
        DIR=`echo ${DIR_DATA}/${DIR_EXPT}/${SUB}/${DATA_TYPE}/${SESS}`
        
        # make eroded white matter mask
        if [ ! -f ${DIR}/anat_wm_ero.nii.gz ]; then
            3dcalc \
                -a ${DIR}/anat_aparc_reg.nii.gz \
                -expr "equals(a,2)  + \
                       equals(a,7)  + \
                       equals(a,41) + \
                       equals(a,46) + \
                       equals(a,251)+ \
                       equals(a,252)+ \
                       equals(a,253)+ \
                       equals(a,254)+ \
                       equals(a,255)" \
                -prefix ${DIR}/anat_wm.nii.gz

            3dcalc \
                -a ${DIR}/anat_wm.nii.gz \
                -b a+i -c a-i -d a+j -e a-j -f a+k -g a-k \
                -expr 'a*(1-amongst(0,b,c,d,e,f,g))' \
                -prefix ${DIR}/anat_wm_ero.nii.gz
        fi

        # make eroded ventricle mask
        if [ ! -f ${DIR}/anat_vent_ero.nii.gz ]; then
            3dcalc \
                -a ${DIR}/anat_aparc_reg.nii.gz \
                -expr 'equals(a,4) + equals(a,43)' \
                -prefix ${DIR}/anat_vent.nii.gz

            3dcalc \
                -a ${DIR}/anat_aparc_reg.nii.gz \
                -expr "equals(a,10) + \
                       equals(a,11) + \
                       equals(a,26) + \
                       equals(a,49) + \
                       equals(a,50) + \
                       equals(a,58)" \
                -prefix ${DIR}/anat_tmp_nonvent.nii.gz

            3dcalc \
                -a ${DIR}/anat_tmp_nonvent.nii.gz \
                -b a+i -c a-i -d a+j -e a-j -f a+k -g a-k \
                -expr 'amongst(1,a,b,c,d,e,f,g)' \
                -prefix ${DIR}/anat_tmp_nonvent_dia.nii.gz

            3dcalc \
                -a ${DIR}/anat_vent.nii.gz \
                -b ${DIR}/anat_tmp_nonvent_dia.nii.gz \
                -expr 'a-step(a*b)' \
                -prefix ${DIR}/anat_vent_ero.nii.gz
        fi

        # make frontal vs temporal mask
        if [ ! -f ${DIR}/anat_temporal_vs_frontal.nii.gz ]; then
            3dcalc \
                -a ${DIR}/anat_aparc_reg.nii.gz \
                -short \
                -expr "1 * equals(a,18) + \
                       2 * equals(a,17) + \
                       5 * equals(a,54) + \
                       6 * equals(a,53) + \
                       3 *(equals(a,1001) + \
                           equals(a,1006) + \
                           equals(a,1007) + \
                           equals(a,1009) + \
                           equals(a,1015) + \
                           equals(a,1016) + \
                           equals(a,1030) + \
                           equals(a,1033) + \
                           equals(a,1034))+ \
                       4 *(equals(a,1002) + \
                           equals(a,1003) + \
                           equals(a,1012) + \
                           equals(a,1014) + \
                           equals(a,1018) + \
                           equals(a,1019) + \
                           equals(a,1020) + \
                           equals(a,1024) + \
                           equals(a,1026) + \
                           equals(a,1027) + \
                           equals(a,1028) + \
                           equals(a,1032))+ \
                       7 *(equals(a,2001) + \
                           equals(a,2006) + \
                           equals(a,2007) + \
                           equals(a,2009) + \
                           equals(a,2015) + \
                           equals(a,2016) + \
                           equals(a,2030) + \
                           equals(a,2033) + \
                           equals(a,2034))+ \
                       8 *(equals(a,2002) + \
                           equals(a,2003) + \
                           equals(a,2012) + \
                           equals(a,2014) + \
                           equals(a,2018) + \
                           equals(a,2019) + \
                           equals(a,2020) + \
                           equals(a,2024) + \
                           equals(a,2026) + \
                           equals(a,2027) + \
                           equals(a,2028) + \
                           equals(a,2032))" \
                -prefix ${DIR}/anat_temporal_vs_frontal.nii.gz
        fi

        # make grey matter mask
        if [ ! -f ${DIR}/anat_gm.nii.gz ]; then
            3dcalc \
                -a ${DIR}/anat_aparc_reg.nii.gz \
                -short \
                -expr 'step(a-1000)*step(1036-a)+step(a-2000)*step(2036-a)' \
                -prefix ${DIR}/anat_gm.nii.gz
        fi

        # make dialated brain mask
        if [ ! -f ${DIR}/anat_EPI_mask_dia.nii.gz ]; then
            3dcalc \
                -a ${DIR}/anat_EPI_mask.nii.gz \
                -b a+i -c a-i -d a+j -e a-j -f a+k -g a-k \
                -expr 'amongst(1,a,b,c,d,e,f,g)' \
                -prefix ${DIR}/anat_EPI_mask_dia.nii.gz
        fi

        # make brainstem mask
        if [ ! -f ${DIR}/anat_bstem.nii.gz ]; then
            3dcalc \
                -a ${DIR}/anat_aparc_reg.nii.gz \
                -expr "equals(a,8)  + \
                       equals(a,47) + \
                       equals(a,16) + \
                       equals(a,12) + \
                       equals(a,13) + \
                       equals(a,26) + \
                       equals(a,51) + \
                       equals(a,52) + \
                       equals(a,17) + \
                       equals(a,18) + \
                       equals(a,53) + \
                       equals(a,54) + \
                       equals(a,58) + \
                       equals(a,28) + \
                       equals(a,60)" \
                -prefix ${DIR}/anat_bstem.nii.gz
        fi

        # make ero draining vessel mask
        if [ ! -f ${DIR}/anat_dv_ero.nii.gz ]; then
            3dcalc \
                -a ${DIR}/anat_EPI_mask.nii.gz \
                -b ${DIR}/anat_gm.nii.gz \
                -c ${DIR}/anat_wm.nii.gz \
                -d ${DIR}/anat_vent.nii.gz \
                -e ${DIR}/anat_tmp_nonvent.nii.gz \
                -f ${DIR}/anat_bstem.nii.gz \
                -expr '(a-b-c-d-e-f)' \
                -prefix ${DIR}/anat_dv.nii.gz
            
            3dcalc \
                -a ${DIR}/anat_dv.nii.gz \
                -b a+i -c a-i -d a+j -e a-j -f a+k -g a-k \
                -expr 'a*(1-amongst(0,b,c,d,e,f,g))' \
                -prefix ${DIR}/anat_dv_ero.nii.gz
        fi

        DIR_RUNS=`ls -d -- ${DIR}/RUN*`
        for RUN in ${DIR_RUNS}; do
            NUM=`basename ${RUN} | sed 's/[^0-9]//g'`

            # create local white matter regressors (+ lagged)
            if [ ! -f ${DIR}/PARAMS/lag.wm_local15.${NUM}.nii.gz ]; then
                3dLocalstat \
                    -prefix ${DIR}/PARAMS/wm_local15.${NUM}.nii.gz \
                    -nbhd 'SPHERE(15)' \
                    -stat mean \
                    -mask ${DIR}/anat_wm_ero.nii.gz \
                    -use_nonmask ${DIR}/func_scaled.${NUM}.nii.gz

                3dTcat \
                    -prefix ${DIR}/PARAMS/lag.wm_local15.${NUM}.nii.gz \
                    ${DIR}/PARAMS/wm_local15.${NUM}.nii.gz'[0]' \
                    ${DIR}/PARAMS/wm_local15.${NUM}.nii.gz'[0..$]'
            fi

            # make ventricle time series (+ lagged)
            if [ ! -f ${DIR}/PARAMS/lag.vent.${NUM}.1D ]; then
                3dmaskave \
                    -q -mask ${DIR}/anat_vent_ero.nii.gz \
                    ${DIR}/func_scaled.${NUM}.nii.gz > \
                    ${DIR}/PARAMS/vent.${NUM}.1D
                
                1dcat \
                    ${DIR}/PARAMS/vent.${NUM}.1D'{0}' > \
                    ${DIR}/PARAMS/lag.vent.${NUM}.1D
                
                1dcat \
                    ${DIR}/PARAMS/vent.${NUM}.1D'{0..$}' >> \
                    ${DIR}/PARAMS/lag.vent.${NUM}.1D
            fi

            # make draining vessel time series (+ lagged)
            if [ ! -f ${DIR}/PARAMS/lag.dv.${NUM}.1D ]; then
                3dmaskave \
                    -q -mask ${DIR}/anat_dv_ero.nii.gz \
                    ${DIR}/func_scaled.${NUM}.nii.gz > \
                    ${DIR}/PARAMS/dv.${NUM}.1D
                
                1dcat \
                    ${DIR}/PARAMS/dv.${NUM}.1D'{0}' > \
                    ${DIR}/PARAMS/lag.dv.${NUM}.1D
                
                1dcat \
                    ${DIR}/PARAMS/dv.${NUM}.1D'{0..$}' >> \
                    ${DIR}/PARAMS/lag.dv.${NUM}.1D
            fi
        done
        rm ${DIR}/anat_tmp*
    done
done
cd ${DIR_PIPE}
