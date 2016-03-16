#!/usr/bin/env python

def run(dir_data, expt, mode):
    output = ''
    print('\nAdding subject-wise motion QC to the outputs.')
    line = '. ${{DIR_PIPE}}/modules/qc/qc_motionind {} {} {} ${{ID}}'.format(dir_data, expt, mode)
    return line, output

