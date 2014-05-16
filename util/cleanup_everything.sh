#!/bin/bash

## Warn the user and ask to continue
read -p "About to delete processed functional data! Press Y/y to continue." \
        -n 1 -r
echo # (optional) move to a new line
if [[ $REPLY =~ ^[Yy]$ ]]
then

## Let's Dance https://www.youtube.com/watch?v=N4d7Wp9kKjA ##
for SUB in ${SUBJECTS}; do

    # remove any multisession files
    rm ${DIR_DATA}/${DIR_EXPT}/${SUB}/${DATA_TYPE}/func_* >& /dev/null

    DIR_SESS=`ls -d -- ${DIR_DATA}/${DIR_EXPT}/${SUB}/${DATA_TYPE}/*/`
    for SESS in ${DIR_SESS}; do
        
        # remove per-session files
        rm -r ${SESS}/PARAMS/ >& /dev/null
        rm ${SESS}/anat_* >& /dev/null
        rm ${SESS}/func_* >& /dev/null
        rm ${SESS}/mat_* >& /dev/null
        rm ${SESS}/reg_* >& /dev/null
        
    done

    DIR_SESS=`ls -d -- ${DIR_DATA}/${DIR_EXPT}/${SUB}/T1/*/`
    for SESS in ${DIR_SESS}; do

        rm ${SESS}/reg_* >& /dev/null
        rm ${SESS}/anat_* >& /dev/null

    done
done
fi
# JDV