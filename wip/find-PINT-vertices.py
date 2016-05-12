#!/usr/bin/env python
"""
Beta version of script find PINT (Personal Instrisic Network Topology)

Usage:
  find-PINT-vertices.py [options] --func <func.dtseries.nii> --surf-L <surf.gii> --surf-R <surf.gii> --input-vertices FILE --outputcsv FILE
Arguments:
    --func <func.dtseries.nii>  Paths to directory source image
    --surf-L <surface.gii>      Path to template for the ROIs of network regions
    --surf-R <surface.gii>      Surface file .surf.gii to read coordinates from
    --input-vertices FILE       Table of template vertices from which to Start
    --outputcsv FILE            Output csv file

Options:
  -v,--verbose             Verbose logging
  --debug                  Debug logging in Erin's very verbose style
  -n,--dry-run             Dry run
  -h,--help                Print help

DETAILS
TBA

Written by Erin W Dickie, April 2016
"""
from epitome.docopt import docopt

import numpy as np
import nibabel as nib
import os
import sys
import tempfile
import shutil
import subprocess
import pandas as pd
import nibabel.gifti.giftiio

arguments     = docopt(__doc__)
func          = arguments['--func']
surfL         = arguments['--surf-L']
surfR         = arguments['--surf-R']
origcsv       = arguments['--input-vertices']
outputfile    = arguments['--outputcsv']
VERBOSE       = arguments['--verbose']
DEBUG         = arguments['--debug']
DRYRUN        = arguments['--dry-run']


#mkdir a tmpdir for the
global RADIUS_SAMPLING
global RADIUS_SEARCH
RADIUS_SAMPLING = 5
RADIUS_SEARCH = 20
tmpdir = tempfile.mkdtemp()

###
if DEBUG: print(arguments)
### Erin's little function for running things in the shell
def docmd(cmdlist):
    "sends a command (inputed as a list) to the shell"
    if DEBUG: print ' '.join(cmdlist)
    if not DRYRUN: subprocess.call(cmdlist)

## measuring distance
def calc_surf_distance(surf, orig_vertex, target_vertex, radius_search, tmpdir=tmpdir):
    '''
    uses wb_command -surface-geodesic-distance command to measure
    distance between two vertices on the surface
    '''
    surf_distance = os.path.join(tmpdir, "distancecalc.shape.gii")
    docmd(['wb_command', '-surface-geodesic-distance',
            surf, str(orig_vertex), surf_distance,
            '-limit', str(radius_search)])
    distances = epi.utilities.load_gii_data(surf_distance)
    distance = distances[target_vertex,0]
    return(distance)

def calc_distance_column(df, orig_vertex_col, target_vertex_col,distance_outcol,
                         radius_search = RADIUS_SEARCH,
                         surfL = surfL, surfR = surfR):
    df.loc[:,distance_outcol] = -99.9
    for idx in df.index.tolist():
        orig_vertex = df.loc[idx, orig_vertex_col]
        target_vertex = df.loc[idx, target_vertex_col]
        hemi = df.loc[idx,'hemi']
        if hemi == "L":
            df.loc[idx, distance_outcol] = calc_surf_distance(surfL, orig_vertex,
                                            target_vertex, radius_search, tmpdir)
        if hemi == "R":
            df.loc[idx, distance_outcol] = calc_surf_distance(surfR, orig_vertex,
                                            target_vertex, radius_search, tmpdir)
    return df

def load_surfaceonly(filename, tempdir = tmpdir):
    '''
    separate a cifti file into surfaces,
    then loads and concatenates the surface data
    '''
    ## separate the cifti file into left and right surfaces
    L_data_surf=os.path.join(tempdir, 'Ldata.func.gii')
    R_data_surf=os.path.join(tempdir, 'Rdata.func.gii')
    docmd(['wb_command','-cifti-separate', filename, 'COLUMN',
        '-metric', 'CORTEX_LEFT', L_data_surf,
        '-metric', 'CORTEX_RIGHT', R_data_surf])

    ## load both surfaces and concatenate them together
    Ldata = epi.utilities.load_gii_data(L_data_surf)
    Rdata = epi.utilities.load_gii_data(R_data_surf)

    return Ldata, Rdata

def roi_surf_data(df, vertex_colname, surf, hemisphere, roi_radius, tmpdir = tmpdir):
    '''
    uses wb_command -surface-geodesic-rois to build rois (3D files)
    then load and collasp that into 1D array
    '''
    ## right the L and R hemisphere vertices from the table out to temptxt
    vertex_list = os.path.join(tmpdir, 'vertex_list.txt')
    df.loc[df.hemi == hemisphere, vertex_colname].to_csv(vertex_list,sep='\n',index=False)

    ## from the temp text build - func masks and target masks
    roi_surf = os.path.join(tmpdir,'roi_surf.func.gii')
    docmd(['wb_command', '-surface-geodesic-rois', surf,
        str(roi_radius),  vertex_list, roi_surf,
        '-overlap-logic', 'EXCLUDE'])

    vlabels = df[df.hemi == hemisphere].roiidx.tolist()
    rois_data = epi.utilities.load_gii_data(roi_surf)
    rois_data = np.multiply(rois_data, vlabels)

    return rois_data

def rois_bilateral(df, vertex_colname, roi_radius, surfL = surfL, surfR = surfR):
    '''
    runs roi_surf_data for both surfaces and combines them to one numpy array
    '''
    rois_L = roi_surf_data(df, vertex_colname, surfL, 'L', roi_radius)
    rois_R = roi_surf_data(df, vertex_colname, surfR, 'R', roi_radius)
    rois = np.vstack((rois_L, rois_R))
    return rois

## loading the dataframe
df = pd.read_csv(origcsv)
df.loc[:,'roiidx'] = pd.Series(np.arange(1,len(df.index)+1), index=df.index)

## load the func data
func_dataL, func_dataR = load_surfaceonly(func, tmpdir)
num_Lverts = func_dataL.shape[0]
func_data = np.vstack((func_dataL, func_dataR))

vertex_incol = 'vertex'
iter_num = 0
max_distance = 10

while iter_num < 10 or max_distance > 1:
    vertex_outcol = 'vertex_{}'.format(iter_num)
    distance_outcol = 'dist_{}'.format(iter_num)
    df.loc[:,vertex_outcol] = -999
    df.loc[:,distance_outcol] = -99.9

    ## load the sampling data
    sampling_rois = rois_bilateral(df, vertex_incol, RADIUS_SAMPLING)

    ## load the search data
    search_rois = rois_bilateral(df, vertex_incol, RADIUS_SEARCH)

    for idx in df.index.tolist():
        vlabel = df.loc[idx,'roiidx']
        network = df.loc[idx,'NETWORK']
        hemi = df.loc[idx,'hemi']
        orig_vertex = df.loc[idx, vertex_incol]
        netlabels = list(set(df[df.NETWORK == network].roiidx.tolist()) - set([vlabel]))
        netseeds = []
        for netlabel in netlabels:
            netseeds = np.hstack((netseeds,np.where(sampling_rois == netlabel)[0]))
        meants = np.mean(func_data[netseeds.astype(int), :], axis=0)
        idx_mask = np.where(search_rois == vlabel)[0]

        # create output array
        seed_corrs = np.zeros(func_data.shape[0])

        # loop through each time series, calculating r
        for i in np.arange(len(idx_mask)):
            seed_corrs[idx_mask[i]] = np.corrcoef(meants, func_data[idx_mask[i], :])[0][1]
        peakvert = np.argmax(seed_corrs, axis=0)
        if hemi =='R': peakvert = peakvert - num_Lverts
        df.loc[idx,vertex_outcol] = peakvert

    ## calc the distances
    df  = calc_distance_column(df, vertex_incol, vertex_outcol, distance_outcol)

    ## print the max distance as things continue..
    max_distance = max(df[distance_outcol])
    print('Iteration {} max distance: {}'.format(iter_num, max_distance))
    vertex_incol = vertex_outcol
    iter_num += 1

## calc a final distance column
df.loc[:,"ivertex"] = df.loc[:,vertex_outcol]
df  = calc_distance_column(df, 'vertex', 'ivertex',
                        'distance', 250)

cols_to_export = c('hemi','NETWORK','roiidx','vertex','ivertex','distance')

df.to_csv(outputfile, columns = cols_to_export, index = False)

#get rid of the tmpdir
shutil.rmtree(tmpdir)
