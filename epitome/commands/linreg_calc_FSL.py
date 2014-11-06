#!/usr/bin/env python

import copy
import epitome.commands as cmd

def linreg_calc_FSL(input_name):
    output = copy.copy(input_name) # return output unharmed

    print('\nCalculating linear registration pathways.')

    try:
        # get the data-quality option
        print('\nSelect data quality:')
        data_quality = ['low', 'high']
        quality = cmd.utils.selector_list(data_quality)

        # set the cost function option
        print('\nRegistration cost function: (see FSL FLIRT for help)')
        cost_fxns = {'mutualinfo' : '= Mutual Information [H(b)+H(s)-H(b,s)]', 
                     'leastsq' : '=  Least Squares [Pearson Correlation]',
                     'corratio' : '= Correlation Ratio', 
                     'normcorr' : '= Normalized Correlation', 
                     'labeldiff' : '= FSL magic!',
                     'bbr' : '= FSL magic!'}

        cost = cmd.utils.selector_dict(cost_fxns)

        # get registration degrees of freedom
        print('\nRegistration degrees of freedom: (see FSL FLIRT for help)')
        degrees_of_freedom = [6, 7, 9, 12]
        reg_dof = cmd.utils.selector_list(degrees_of_freedom)

    # if we messed any of these up, we return None
    except ValueError as ve:
        return '', None

    # otherwise we print the command and return it
    line = ('. ${DIR_PIPE}/epitome/modules/pre/linreg_calc_FSL ' +
                                               str(quality) + ' ' +
                                               str(cost) + ' ' +
                                               str(reg_dof))
    return line, output