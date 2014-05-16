#!/usr/bin/env python

import os
import sys
import datetime

import nibabel as nib
import numpy as np
from scipy import ndimage as nd

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import ypp_utilities

def reg_check(path, expt, mode):
    """
    Prints the central slice of the T1 and co-registered + deskulled EPI, 
    including an edge-detected version of the T1 (requires AFNI).
    """

    # get subject numbers
    subjects = ypp_utilities.get_subj(os.path.join(path, expt))

    # loop through all subjects
    pdf = PdfPages(os.path.join(path, expt, 'qc_reg_EPI_to_T1.pdf'))
    for subj in subjects:

        edge = os.path.join(path, expt, subj, 'T1/SESS01/anat_T1_edge.nii.gz')
        anat = os.path.join(path, expt, subj, 'T1/SESS01/anat_T1_brain.nii.gz')
        reg = os.path.join(path, expt, subj, mode, 'SESS01/reg_EPI_to_T1.nii.gz')

        # create edge dataset if it doesn't exist
        print 'working on subject ' + str(subj)
        if os.path.isfile(edge) == False:
            os.system('3dedge3 -input ' + anat + ' -prefix ' + edge)

        # load in data
        edge = nib.load(edge).get_data()
        anat = nib.load(anat).get_data()
        reg = nib.load(reg).get_data()

        # reorient the data to radiological
        edge = np.transpose(edge, (2,0,1))
        edge = np.rot90(edge, 2)
        anat = np.transpose(anat, (2,0,1))
        anat = np.rot90(anat, 2)
        reg = np.transpose(reg, (2,0,1))
        reg = np.rot90(reg, 2)

        # get size ratio between over + underlay
        dsfactor = [a/float(r) for a,r in zip(anat.shape, reg.shape)]
        # match over + underlay dimensions
        reg_to_anat = nd.interpolation.zoom(reg, zoom=dsfactor)
        # set small values in overlay to be transparent
        reg_to_anat = np.ma.masked_where(reg_to_anat < 50, reg_to_anat)
        cmap = plt.cm.Reds
        cmap.set_bad('g', 0)

        # generate the overlay image
        plt.subplot(2,3,1)
        mid = np.round(anat.shape[0] / 2)
        plt.imshow(anat[mid, :, :], cmap=plt.cm.gray,
                                    interpolation='nearest')
        plt.imshow(reg_to_anat[mid, :, :], cmap=cmap,
                                           interpolation='nearest',
                                           alpha=0.5)
        plt.axis('off')

        plt.subplot(2,3,2)
        mid = np.round(anat.shape[1] / 2)
        plt.imshow(anat[:, mid, :], cmap=plt.cm.gray,
                                    interpolation='nearest')
        plt.imshow(reg_to_anat[:, mid, :], cmap=cmap,
                                           interpolation='nearest',
                                           alpha=0.5)
        plt.axis('off')

        plt.subplot(2,3,3)
        mid = np.round(anat.shape[2] / 2)
        plt.imshow(anat[:, :, mid], cmap=plt.cm.gray,
                                    interpolation='nearest')
        plt.imshow(reg_to_anat[:, :, mid], cmap=cmap,
                                           interpolation='nearest',
                                           alpha=0.5)
        plt.axis('off')

        # set zeros in edge to be transparent
        edge = np.ma.masked_where(edge == 0, edge)
        cmap = plt.cm.winter
        cmap.set_bad('g', 0)

        # generate the edge image
        plt.subplot(2,3,4)
        mid = np.round(edge.shape[0] / 2)
        plt.imshow(reg_to_anat[mid, :, :], cmap=plt.cm.gray,
                                           interpolation='nearest')
        plt.imshow(edge[mid, :, :], cmap=cmap,
                                    interpolation='nearest')
        plt.axis('off')

        plt.subplot(2,3,5)
        mid = np.round(edge.shape[1] / 2)
        plt.imshow(reg_to_anat[:, mid, :], cmap=plt.cm.gray,
                                           interpolation='nearest')
        plt.imshow(edge[:, mid, :], cmap=cmap,
                                    interpolation='nearest')
        plt.axis('off')

        plt.subplot(2,3,6)
        mid = np.round(edge.shape[2] / 2)
        plt.imshow(reg_to_anat[:, :, mid], cmap=plt.cm.gray,
                                           interpolation='nearest')
        plt.imshow(edge[:, :, mid], cmap=cmap,
                                    interpolation='nearest')
        plt.axis('off')

        plt.suptitle(str(expt) + ' ' + str(mode) + ': ' + str(subj))
        plt.tight_layout()
        plt.savefig(pdf, format='pdf')
        plt.close()

    # Add some metadata and close the PDF object
    d = pdf.infodict()
    d['Title'] = 'Quality Control: Registration of the EPI template to the T1'
    d['Author'] = u'Joseph D Viviano\xe4nen'
    d['Subject'] = 'Quality Control'
    d['Keywords'] = 'QC registration EPI T1'
    d['CreationDate'] = datetime.datetime.today()
    d['ModDate'] = datetime.datetime.today()
    pdf.close()

if __name__ == "__main__":
	reg_check(sys.argv[1], sys.argv[2], sys.argv[3])
