#!/usr/bin/env python

def run(dir_data, expt, mode):
    output = ''
    print('\nAdding NIFTI dimension-checking QC to the outputs.')
    line = '. ${DIR_PIPE}/epitome/modules/qc/qc_runs {} {}'.format(dir_data, expt)
    return line, output

