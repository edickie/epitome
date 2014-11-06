#!/usr/bin/env python

import epitome.commands as cmd

def nonlinreg_EPI2MNI_AFNI(input_name):
    output = 'MNI'

    # give us some feedback
    print('\nNonlinearly re-sampling input EPI data to MNI space using AFNI.')

    try:
        # get the reslice dimensions
        print('\nSelect target dimensions (isotropic mm):')
        dims = cmd.utils.selector_float()

    # if we messed any of these up, we return None
    except ValueError as ve:
        return '', None

    # otherwise we print the command and return it
    line = ('. ${DIR_PIPE}/epitome/modules/pre/linreg_EPI2MNI_FSL ' +
                                              str(input_name) + ' ' +
                                              str(dims))
    return line, output