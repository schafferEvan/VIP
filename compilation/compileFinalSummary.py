#!/usr/bin/env python

""" aggregate data from processing pipeline and export final output to npz files to share

@author evan schaffer
"""
# Created July 17 2019


import os
import sys
import glob
import numpy as np
from scipy import sparse, io


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
    matPost=io.loadmat(expDir+'/post_fromYcc.mat') 
    matRaw=io.loadmat(expDir+'/F.mat') 
else:
    matPost=io.loadmat(expDir+'/post_fromRcc.mat') 
    matRaw=io.loadmat(expDir+'/F_fromRed.mat') 

matSum = io.loadmat(expDir+'/Ysum.mat') 
matbeh = io.loadmat(expDir+'alignedBehavAndStim.mat')

trialFlag=matRaw['trialFlag']
dFF=matPost['dOO']
A=matRaw['A']

Ysum=matSum['Ysum']
dims=np.shape(Ysum)

time = matbeh['time']
ball = matbeh['alignedBehavior']['legVar'][0][0] 
dlc = matbeh['alignedBehavior']['dlcData'][0][0] 
goodIds = np.flatnonzero( matPost['goodIds'] )

dFF = matPost['dOO']
dFF = dFF[goodIds,:]
A = A[:,goodIds]

# compress Ysum
M = np.max(Ysum)
m = np.min(Ysum)
im = np.uint16(2**16*(Ysum-m)/(M-m))




# save stuff --------------------------------------

expNameHandle=expID.replace('/','_')
saveHandle = savePath+expNameHandle

np.savez( saveHandle+'.npz', time=time, trialFlag=trialFlag, dFF=dFF, ball=ball, dlc=dlc, dims=dims) #dims=dims, A=A, 
sparse.save_npz(saveHandle+'_A.npz', A)

np.savez( saveHandle+'_ref_im.npz', im=im) #dims=dims, A=A, 
