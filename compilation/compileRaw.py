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
    fromGreenCC = False    # use ROIs parsed from green image (false -> use red) # DEPRECATED. ALWAYS DO BOTH NOW


# load stuff --------------------------------------
# matRaw_fromGreen=io.loadmat(expDir+'/F.mat') 
matRaw_fromRed=io.loadmat(expDir+'/F_fromRed.mat') 

matSum = io.loadmat(expDir+'/Ysum.mat') 

Y=matRaw_fromRed['F']
R=matRaw_fromRed['FR']
# YfY=matRaw_fromGreen['F']

trialFlag=matRaw_fromRed['trialFlag']
A=matRaw_fromRed['A']
# AfY=matRaw_fromGreen['A']

Ysum=matSum['Ysum']
dims=np.shape(Ysum)

# get regionProps
matcc = io.loadmat(expDir+'/cc.mat')
try:
    rprops = matcc['regionProps']['blobStats'][0][0]
    centroids = np.zeros(rprops.shape[0],3)
    for i in range(rprops.shape[0]):
        centroids[i,:] = rprops[i]['Centroid'][0]
except:
    print('Warning: No centroids from Red')
    centroids = []

try:
    rprops_g = matcc['regionProps_green']['blobStats'][0][0]
    centroids_fromGreen = np.zeros(rprops_g.shape[0],3)
    for i in range(rprops_g.shape[0]):
        centroids_fromGreen[i,:] = rprops_g[i]['Centroid'][0]
except:
    print('Warning: No centroids from Green')
    centroids_fromGreen = []


if os.path.isfile(expDir+'alignedBehavAndStim.mat'):           
    try:
        matbeh = io.loadmat(expDir+'alignedBehavAndStim.mat')
        time = matbeh['time'].T
        ball = matbeh['alignedBehavior']['legVar'][0][0].T
        dlc = matbeh['alignedBehavior']['dlcData'][0][0].T
        stim = matbeh['alignedBehavior']['stim'][0][0].T
        drink = matbeh['alignedBehavior']['drink'][0][0].T
        print('v7')
    except:
        print('v7.3')
        matbeh = h5py.File(expDir+'alignedBehavAndStim.mat')
        time = matbeh['time']
        ball = matbeh['alignedBehavior']['legVar']
        dlc = matbeh['alignedBehavior']['dlcData']
        stim = matbeh['alignedBehavior']['stim']
        drink = matbeh['alignedBehavior']['drink']
else:
    print('no aligned behavior found')
    time = []
    ball = []
    dlc = []
    stim = []
    drink = []         

# check if behavior label file exists
try:
    label_file = h5py.File(expDir + 'alignedLabels.mat')
    states = label_file['labelsAligned']
    print('found behavior labels')
except:
    states = []
    print('no behavior labels')

# check if PID timing file exists
try:
    PID_file = h5py.File(expDir + 'alignedPID.mat')
    PIDAligned = PID_file['PIDAligned']
    print('found aligned PID trace')
except:
    PIDAligned = []
    print('no PID trace found')    

# compress Ysum
M = np.max(Ysum)
m = np.min(Ysum)
im = np.uint16(2**16*(Ysum-m)/(M-m))


# get scanrate
infodir = expDir+'info/'
infols=glob.glob(infodir+'f*.mat')


info = io.loadmat(infols[0])
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

np.savez( saveHandle+'_raw.npz', scanRate=scanRate, time=time, trialFlag=trialFlag, Y=Y, R=R, ball=ball, dlc=dlc, stim=stim, drink=drink,
            dims=dims, im=im, tot_um_x=tot_um_x, tot_um_y=tot_um_y, tot_um_z=tot_um_z, states=states, PIDAligned=PIDAligned, 
            centroids=centroids, centroids_fromGreen=centroids_fromGreen) 
sparse.save_npz(saveHandle+'_A_raw.npz', A)

# io.savemat(saveHandle+'_raw.mat',{
#     'time':time, 'trialFlag':trialFlag, 'Y':Y, 'R':R, 'ball':ball, 'dlc':dlc, 'dims':dims, 'A':A, 'im':im})

# f = h5py.File(saveHandle+"_raw.hdf5", "w")
# f.create_dataset('Y', Y)
# f.create_dataset('R', R)