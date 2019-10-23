#!/usr/bin/env python

""" aggregate data from processing pipeline

@author evan schaffer
"""
# Created July 17 2019


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
    fromGreenCC = mynargs[4]=='True'

else:
    expID = '2019_06_26_Nsyb_NLS6s_walk/fly2'
    expDir = '/Volumes/SCAPEdata1/finalData/'+expID+'/Yproj/'
    savePath = expDir #'/Users/evan/Downloads/' #'/Users/evan/Dropbox/_AxelLab/_flygenvectors_dataShare/'
    fromGreenCC = False    # use ROIs parsed from green image (false -> use red)


# load stuff --------------------------------------

if fromGreenCC:
    matRaw=io.loadmat(expDir+'/F.mat') 
else:
    matRaw=io.loadmat(expDir+'/F_fromRed.mat') 

matSum = io.loadmat(expDir+'/Ysum.mat') 

Y=matRaw['F']
R=matRaw['FR']

trialFlag=matRaw['trialFlag']
A=matRaw['A']

Ysum=matSum['Ysum']
dims=np.shape(Ysum)

# # get regionProps
# matcc = io.loadmat(expDir+'/cc.mat')
# try:
#     rprops = matcc['regionProps']['blobStats'][0][0]
#     centroids = np.zeros(rprops.shape[0],3)
#     for i in range(rprops.shape[0]):
#         centroids[i,:] = rprops[i]['Centroid'][0]
# except:
#     centroids = []


try:
    matbeh = io.loadmat(expDir+'alignedBehavAndStim.mat')
    time = matbeh['time'].T
    ball = matbeh['alignedBehavior']['legVar'][0][0].T
    dlc = matbeh['alignedBehavior']['dlcData'][0][0].T
    print('v7')
except:
    print('v7.3')
    matbeh = h5py.File(expDir+'alignedBehavAndStim.mat')
    time = matbeh['time']
    ball = matbeh['alignedBehavior']['legVar']
    dlc = matbeh['alignedBehavior']['dlcData']

# check if behavior label file exists
try:
    label_file = h5py.File(expDir + 'alignedLabels.mat')
    states = label_file['labelsAligned']
    print('found behavior labels')
except:
    states = []
    print('no behavior labels')

# check if registered pointset file exists
try:
    point_file = io.loadmat(expDir + 'registered_pointset.mat')
    raw = point_file['raw']
    aligned = point_file['aligned']
    print('found registered pointset')
except:
    try:
        point_file = h5py.File(expDir + 'registered_pointset.mat')
        raw = point_file['raw']
        aligned = point_file['aligned']
        print('found registered pointset')
    except:
        raw = []
        aligned = []
        print('no registered pointset')

# compress Ysum
M = np.max(Ysum)
m = np.min(Ysum)
im = np.uint16(2**16*(Ysum-m)/(M-m))


# get scanrate
infodir = expDir+'info/'
infols=os.listdir(infodir)
info = io.loadmat(infodir+infols[0])
scanRate = info['info']['daq'][0][0]['scanRate'][0][0][0][0]
#pdb.set_trace()


# get true dimensions of image in um
try:
    x_umPerPix = info['info']['GUIcalFactors'][0][0]['x_umPerPix'][0][0][0][0]
    tot_um_x = dims[0]*x_umPerPix
except:
    xK_umPerVolt = info['info']['GUIcalFactors'][0][0]['xK_umPerVolt'][0][0][0][0]
    scanAngle = info['info']['daq'][0][0]['scanAngle'][0][0][0][0]
    pixelsPerLine = info['info']['daq'][0][0]['pixelsPerLine'][0][0][0][0]
    tot_um_x = dims[0]*xK_umPerVolt*scanAngle/(pixelsPerLine-1);

y_umPerPix = info['info']['GUIcalFactors'][0][0]['y_umPerPix'][0][0][0][0]
tot_um_y = dims[1]*y_umPerPix
z_umPerPix = info['info']['GUIcalFactors'][0][0]['z_umPerPix'][0][0][0][0]
tot_um_z = dims[2]*z_umPerPix

# save stuff --------------------------------------

expNameHandle=expID.replace('/','_')
saveHandle = savePath+expNameHandle

np.savez( saveHandle+'_raw.npz', scanRate=scanRate, time=time, trialFlag=trialFlag, Y=Y, R=R, ball=ball, dlc=dlc, 
            dims=dims, im=im, tot_um_x=tot_um_x, tot_um_y=tot_um_y, tot_um_z=tot_um_z, states=states, 
            original_centroids=raw, aligned_centroids=aligned) 
sparse.save_npz(saveHandle+'_A_raw.npz', A)

# io.savemat(saveHandle+'_raw.mat',{
#     'time':time, 'trialFlag':trialFlag, 'Y':Y, 'R':R, 'ball':ball, 'dlc':dlc, 'dims':dims, 'A':A, 'im':im})

# f = h5py.File(saveHandle+"_raw.hdf5", "w")
# f.create_dataset('Y', Y)
# f.create_dataset('R', R)