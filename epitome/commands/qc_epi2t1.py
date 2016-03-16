#!/usr/bin/env python

def run(dir_data, expt, mode):
    output = ''
    print('\nAdding EPI-to-T1 registration checking QC to the outputs.')
    line = '. ${{DIR_PIPE}}/epitome/modules/qc/qc_epi2t1 {} {} {}'.format(dir_data, expt, mode)
    return line, output

