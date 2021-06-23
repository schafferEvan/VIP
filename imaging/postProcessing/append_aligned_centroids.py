#!/usr/bin/env python

""" aggregate data from processing pipeline

@author evan schaffer
"""
# Created Oct 23 2019


import os
import sys
import glob
import numpy as np
from scipy import sparse, io
import h5py
import pdb


mynargs = sys.argv
if (len(mynargs)>1):
    expDir = mynargs[1]
    expID = mynargs[2]
    savePath = mynargs[3]

else:
    expID = '2019_06_26_Nsyb_NLS6s_walk/fly2'
    expDir = '/Volumes/SCAPEdata1/finalData/'+expID+'/Yproj/'
    savePath = expDir #'/Users/evan/Downloads/' #'/Users/evan/Dropbox/_AxelLab/_flygenvectors_dataShare/'

expNameHandle=expID.replace('/','_')
fileHandle = savePath+expNameHandle+'.npz'

f = np.load( fileHandle , allow_pickle=True)
try:
    ptf = io.loadmat( savePath+'registered_pointset.mat')
    print('v7')
except:
    ptf = h5py.File(savePath+'registered_pointset.mat')
    print('v7.3')
# pdb.set_trace()

# sparse red datasets (registered on pan green) have flag for green/red centroid ID
if 'is_red' in ptf:
    centroid_is_red = ptf['is_red']
else:
    centroid_is_red = []

# log volume dims to match template (see make_template.m)
range_buffer = 1.02
template_file = io.loadmat( '/Volumes/SCAPEdata1/finalData/_templates/dorsopost_template.mat' )
dims_in_um_adj = range_buffer*template_file['ref'].max(axis=0)[[1,0,2]]

np.savez( fileHandle, time=f['time'], trialFlag=f['trialFlag'],
        dFF=f['dFF'], dYY=f['dYY'], dRR=f['dRR'], ball=f['ball'], dlc=f['dlc'], 
        beh_labels=f['beh_labels'], stim=f['stim'], drink=f['drink'],
        dims=f['dims'], dims_in_um=dims_in_um_adj, dims_in_um_orig=f['dims_in_um'], im=f['im'], 
        scanRate=f['scanRate'], redTh=f['redTh'], grnTh=f['grnTh'],
        goodIds=f['goodIds'], oIsGood=f['oIsGood'], PIDdata=f['PIDdata'],
        aligned_centroids=ptf['aligned'], centroid_is_red=centroid_is_red ) 
print('--- COMPLETED APPENDING ALIGNED CENTROIDS ---')


