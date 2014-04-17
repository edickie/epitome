#!/bin/bash

## Filter Module: Regresses computed nusiance variables from each run
# Computes detrended nusiance time series
# Fits each run with all nusiance variables
# Computed temporal SNR
# Subtracts noise time series from each voxel, retaining the mean

cd /tmp
for SUB in ${SUBJECTS}; do
    DIR_SESS=`ls -d -- ${DIR_DATA}/${DIR_EXPT}/${SUB}/${DATA_TYPE}/*/`
    for SESS in ${DIR_SESS}; do
        # this used to be here: `ls -d -- ${SESS}/run*/ | sort -n -t n -k 4`
        DIR_RUNS=`ls -d -- ${SESS}/RUN*`
        for RUN in ${DIR_RUNS}; do
            NUM=`basename ${RUN} | sed 's/[^0-9]//g'`

            # detrend datasets and .1D files
            #3dDetrend -prefix - \
            #          -DAFNI_1D_TRANOUT=YES \
            #          -polort ${POLORT} \
            #         # RVT REGRESSOR GOES HERE IN CONDITIONAL LOOP
            if [ ! -f ${SESS}/PARAMS/det.motion.${NUM}.1D ]; then
                3dDetrend \
                    -prefix - \
                    -DAFNI_1D_TRANOUT=YES \
                    -polort ${POLORT} \
                    ${SESS}/PARAMS/motion.${NUM}.1D\' > \
                    ${SESS}/PARAMS/det.motion.${NUM}.1D
            fi
            
            if [ ! -f ${SESS}/PARAMS/det.vent.${NUM}.1D ]; then
                3dDetrend \
                    -prefix - \
                    -DAFNI_1D_TRANOUT=YES \
                    -polort ${POLORT} \
                    ${SESS}/PARAMS/vent.${NUM}.1D\' > \
                    ${SESS}/PARAMS/det.vent.${NUM}.1D
            fi
            
            if [ ! -f ${SESS}/PARAMS/det.lag.vent.${NUM}.1D ]; then
                3dDetrend \
                    -prefix - \
                    -DAFNI_1D_TRANOUT=YES \
                    -polort ${POLORT} \
                    ${SESS}/PARAMS/lag.vent.${NUM}.1D\' > \
                    ${SESS}/PARAMS/det.lag.vent.${NUM}.1D
            fi

            if [ ! -f ${SESS}/PARAMS/det.dv.${NUM}.1D ]; then
                3dDetrend \
                    -prefix - \
                    -DAFNI_1D_TRANOUT=YES \
                    -polort ${POLORT} \
                    ${SESS}/PARAMS/dv.${NUM}.1D\' > \
                    ${SESS}/PARAMS/det.dv.${NUM}.1D
            fi
            
            if [ ! -f ${SESS}/PARAMS/det.lag.dv.${NUM}.1D ]; then
                3dDetrend \
                    -prefix - \
                    -DAFNI_1D_TRANOUT=YES \
                    -polort ${POLORT} \
                    ${SESS}/PARAMS/lag.dv.${NUM}.1D\' > \
                    ${SESS}/PARAMS/det.lag.dv.${NUM}.1D
            fi
            
            if [ ! -f ${SESS}/PARAMS/det.wm_local15.${NUM}.nii.gz ]; then
                3dDetrend \
                    -prefix ${SESS}/PARAMS/det.wm_local15.${NUM}.nii.gz \
                    -polort ${POLORT} \
                    ${SESS}/PARAMS/wm_local15.${NUM}.nii.gz
            fi

            if [ ! -f ${SESS}/PARAMS/det.lag.wm_local15.${NUM}.nii.gz ]; then
                3dDetrend \
                    -prefix ${SESS}/PARAMS/det.lag.wm_local15.${NUM}.nii.gz \
                    -polort ${POLORT} \
                    ${SESS}/PARAMS/lag.wm_local15.${NUM}.nii.gz
            fi

            if [ ! -f ${SESS}/PARAMS/det.global_mean.${NUM}.1D ]; then
                3dDetrend \
                    -prefix - \
                    -DAFNI_1D_TRANOUT=YES \
                    -polort ${POLORT} \
                    ${SESS}/PARAMS/global_mean.${NUM}.1D\' > \
                    ${SESS}/PARAMS/det.global_mean.${NUM}.1D
            fi

            # fits each run with all nusiance variables
            # need if confitional for RVT...
            if [ ! -f ${SESS}/func_noise.${NUM}.nii.gz ]; then
                3dTfitter \
                    -prefix ${SESS}/func_noise_betas.${NUM}.nii.gz \
                    -fitts ${SESS}/func_noise.${NUM}.nii.gz \
                    -polort ${POLORT} \
                    -RHS ${SESS}/func_scaled.${NUM}.nii.gz \
                    -LHS ${SESS}/PARAMS/det.motion.${NUM}.1D \
                         ${SESS}/PARAMS/det.vent.${NUM}.1D \
                         ${SESS}/PARAMS/det.lag.vent.${NUM}.1D \
                         ${SESS}/PARAMS/det.dv.${NUM}.1D \
                         ${SESS}/PARAMS/det.lag.dv.${NUM}.1D \
                         ${SESS}/PARAMS/det.wm_local15.${NUM}.nii.gz \
                         ${SESS}/PARAMS/det.lag.wm_local15.${NUM}.nii.gz \
                         ${SESS}/PARAMS/det.global_mean.${NUM}.1D
            fi

            # compute mean, standard deviation
            if [ ! -f ${SESS}/func_filtered.${NUM}.nii.gz ]; then
                3dTstat \
                    -prefix ${SESS}/func_tmp_mean.${NUM}.nii.gz \
                    -mean ${SESS}/func_scaled.${NUM}.nii.gz
               
                3dTstat \
                    -prefix ${SESS}/func_tmp_stdev.${NUM}.nii.gz \
                    -stdev ${SESS}/func_scaled.${NUM}.nii.gz
                
                # compute temporal SNR
                3dcalc \
                    -a ${SESS}/func_tmp_mean.${NUM}.nii.gz \
                    -b ${SESS}/func_tmp_stdev.${NUM}.nii.gz \
                    -expr 'a/b' \
                    -prefix ${SESS}/func_tSNR.${NUM}.nii.gz
                
                # subtracts nusiances from inputs, retaining the mean
                3dcalc \
                    -float \
                    -a ${SESS}/func_scaled.${NUM}.nii.gz \
                    -b ${SESS}/func_noise.${NUM}.nii.gz \
                    -c ${SESS}/func_tmp_mean.${NUM}.nii.gz \
                    -expr 'a-b+c' \
                    -prefix ${SESS}/func_filtered.${NUM}.nii.gz
            fi
        done
        rm ${SESS}/func_tmp_*
    done
done
cd ${DIR_PIPE}
