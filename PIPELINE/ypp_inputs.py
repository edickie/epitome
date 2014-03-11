#!/usr/bin/env python

# defines paths, subjects, etc.

def init():
    path = '/srv/MRI/WORKING'                    # main directory
    mode = 'REST'
    core = 7

    #expt = 'BEBASD'
    #subj = ['210']

    expt = 'RSFC1'
    subj = ['0023', '0028', '0031', '0035', '0038', 
            '0041', '0044', '0048', '0052', '0058',
            '0062', '0066', '0070', '0073', '0076',
            '0079', '0083', '0024', '0029', '0033',
            '0036', '0039', '0042', '0046', '0050',
            '0054', '0060', '0064', '0067', '0071',
            '0074', '0077', '0081', '0084', '0025',
            '0030', '0034', '0037', '0040', '0043',
            '0047', '0051', '0055', '0061', '0065',
            '0068', '0072', '0075', '0078', '0082']


# expt = 'atol'                                     # experiment directory
# subj = ['301', '302', '304', '305', '307', '308', 
#         '309', '311', '312', '313', '316', '317', 
#         '319', '320', '321', '322', '323', '324',
#         '401', '403', '404', '405', '407', '408', 
#         '409', '411', '412', '413', '414', '415', 
#         '416', '418', '419', '420', '421', '422']

#    expt = 'TRSE'                                     # experiment directory
#    subj = ['1101', '1103', '1104', '1202', '1205', 
#            '1208', '1209', '1210', '1212', '1214', 
#            '1220', '1223', '1306', '1307', '1309', 
#            '1310', '1311', '1313', '1314', '1318', 
#            '1325', '1326', '1328', '1329', '1331', 
#            '1332', '1333', '1336', '1337', '1338', 
#            '1339', '1340', '1343', '1344', '1346', 
#            '1347', '1349', '1350']


# expt = 'trse_enhance'                              # experiment directory
# subj = ['700', '701', '702', '703', '704', '705',   # input participants
#         '706', '707', '708', '709', '710', '712']


    return path, expt, subj, mode, core
