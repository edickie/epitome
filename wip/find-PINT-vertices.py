#!/usr/bin/env python
"""
Beta version of script find PINT (Personal Instrisic Network Topology)

Usage:
  find-PINT-vertices.py [options] --func <func.dtseries.nii> --surf-L <surf.gii> --surf-R <surf.gii> --input-vertices FILE
Arguments:
    --func <func.dtseries.nii>  Paths to directory source image
    --surf-L <surface.gii>      Path to template for the ROIs of network regions
    --surf-R <surface.gii>      Surface file .surf.gii to read coordinates from
    --input-vertices FILE       Table of template vertices from which to Start

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
VERBOSE       = arguments['--verbose']
DEBUG         = arguments['--debug']
DRYRUN        = arguments['--dry-run']

###
if DEBUG: print(arguments)
### Erin's little function for running things in the shell
def docmd(cmdlist):
    "sends a command (inputed as a list) to the shell"
    if DEBUG: print ' '.join(cmdlist)
    if not DRYRUN: subprocess.call(cmdlist)

## measuring distance
def calc_surf_distance(surf, orig_vertex, target_vertex, radius_search, tmpdir):
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

def load_surfaceonly(filename, tempdir):
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

def roi_surf_data(df, surf, hemisphere, roi_radius, tmpdir):
    '''
    uses wb_command -surface-geodesic-rois to build rois (3D files)
    then load and collasp that into 1D array
    '''
    ## right the L and R hemisphere vertices from the table out to temptxt
    vertex_list = os.path.join(tmpdir, 'vertex_list.txt')
    df[df.hemisphere == hemisphere].vertex.to_csv(vertex_list,sep='\n',index=False)

    ## from the temp text build - func masks and target masks
    roi_surf = os.path.join(tmpdir,'roi_surf.func.gii')
    docmd(['wb_command', '-surface-geodesic-rois', surf, vertex_list,
        str(roi_radius), -overlap-logic, 'EXCLUDE',roi_surf])

    vlabels = df[df.hemisphere == hemisphere].roiidx.tolist()
    rois_data = epi.utilities.load_gii_data(roi_surf)
    rois_data = np.multiply(rois_data, vlables)

    return rois_data


#mkdir a tmpdir for the
tmpdir = tempfile.mkdtemp()
radius_sampling='5'
radius_search='10'

## loading the dataframe
df = pd.read_csv(origcsv)
df.loc[:,'roiidx'] = pd.Series(np.arange(1,len(df.index)+1), index=df.index)
df.loc[:,'ivertex'] = -999
df.loc[:,'distance'] = -99.9

## load the func data
func_dataL, func_dataR = load_surfaceonly(func, tmpdir)
num_Lverts = func_dataL.shape[0]
func_data = np.vstack((func_dataL, func_dataR))

## load the sampling data
sampling_rois_L = roi_surf_data(df, surfL, 'L', radius_sampling, tmpdir)
sampling_rois_R = roi_surf_data(df, surfR, 'R', radius_sampling, tmpdir)
sampling_rois = np.vstack((sampling_rois_L, sampling_rois_R))
del sampling_rois_R sampling_rois_L

## load the search data
search_rois_L = roi_surf_data(df, surfL, 'L', radius_search, tmpdir)
search_rois_R = roi_surf_data(df, surfR, 'R', radius_search, tmpdir)
search_rois = np.vstack((search_rois_L, search_rois_R))
del search_rois_L search_rois_R

for idx in df.index.tolist():
    vlabel = df.loc[idx,'roiidx']
    network = df.loc[idx,'NETWORK']
    hemi = df.loc[idx,'hemisphere']
    netlabels = list(set(df[df.NETWORK == network].roiidx.tolist()) - set([vlabel]))
    netseeds = []
    for netlabel in netlabels:
        netseeds = np.hstack(netseeds,np.where(sampling_rois == netlabel)[0])
    meants = np.mean(func_data[netseeds, :], axis=0)
    idx_mask = np.where(search_rois == vlabel)[0]
    # create output array
    seed_corrs = np.zeros(func_data.shape)
    # look through each time series, calculating r
    for i in np.arange(len(idx_mask)):
        seed_corrs[idx_mask[i]] = np.corrcoef(meants, func[idx_mask[i], :])[0][1]
    peakvert = numpy.argmax(seed_corrs, axis=0, out=None)
    if hemi =='R': peakvert = peakvert - num_Lverts
    df.loc[idx,'ivertex'] = peakvert
    if hemi == "L":
        df.loc[idx, 'distance'] = calc_surf_distance(surfL, orig_vertex,
                                        peakvert, radius_search, tmpdir)
    if hemi == "R":
        df.loc[idx, 'distance'] = calc_surf_distance(surfR, orig_vertex,
                                        peakvert, radius_search, tmpdir)                                    

#get rid of the tmpdir
shutil.rmtree(tmpdir)
