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
# import pdb


mynargs = sys.argv
if (len(mynargs)>1):
    expDir = mynargs[1]
    expID = mynargs[2]
    savePath = mynargs[3]
    fromGreenCC = mynargs[4]=='True'

else:
    expDir = '/Volumes/SCAPEdata1/finalData/2019_07_01_Nsyb_NLS6s_walk/fly2/Yproj/'
    expID = '2019_07_01_Nsyb_NLS6s_walk/fly2'
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

# pdb.set_trace()
try:
    matbeh = io.loadmat(expDir+'alignedBehavAndStim.mat')
    time = matbeh['time']
    ball = matbeh['alignedBehavior']['legVar'][0][0] 
    dlc = matbeh['alignedBehavior']['dlcData'][0][0] 
except:
    # ver 7.3
    matbeh = h5py.File(expDir+'alignedBehavAndStim.mat')
    time = matbeh['time']
    ball = matbeh['alignedBehavior']['legVar']
    dlc = matbeh['alignedBehavior']['dlcData']

# compress Ysum
M = np.max(Ysum)
m = np.min(Ysum)
im = np.uint16(2**16*(Ysum-m)/(M-m))




# save stuff --------------------------------------

expNameHandle=expID.replace('/','_')
saveHandle = savePath+expNameHandle

np.savez( saveHandle+'_raw.npz', time=time, trialFlag=trialFlag, Y=Y, R=R, ball=ball, dlc=dlc, dims=dims, im=im) #dims=dims, A=A, 
sparse.save_npz(saveHandle+'_A_raw.npz', A)

# io.savemat(saveHandle+'_raw.mat',{
#     'time':time, 'trialFlag':trialFlag, 'Y':Y, 'R':R, 'ball':ball, 'dlc':dlc, 'dims':dims, 'A':A, 'im':im})

# f = h5py.File(saveHandle+"_raw.hdf5", "w")
# f.create_dataset('Y', Y)
# f.create_dataset('R', R)