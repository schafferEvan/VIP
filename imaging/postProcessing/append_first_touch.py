#!/usr/bin/env python

""" aggregate data from processing pipeline

@author evan schaffer
"""
# Created Sept 10 2020


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
    ptf = io.loadmat( savePath+'sucrose_touch_aligned.mat')
    print('v7')
except:
    ptf = h5py.File(savePath+'sucrose_touch_aligned.mat')
    print('v7.3')
# pdb.set_trace()

s = f['stim']
t = f['time']

# find framestamp of first stim (touch timestamp is stored relative to this point)
s1 = next(x for x, val in enumerate(s) if val > 0) 

# frame stamp of touch is time closest to time of first stim minus touch period
v = np.argmin(abs( t - (t[s1]+ptf['touch_timestamp']) ))

# sucrose_touch is binary vector correctly aligned to other beh, 0 until touch, 1 thereafter
sucrose_touch = np.zeros(s.shape)
sucrose_touch[v:] = 1

# plt.figure(figsize=(16,3))
# plt.plot(t,s)
# plt.plot(t[s1],s[s1],'r.')
# plt.plot(t[v],s[v],'c.')
# plt.plot(t,sucrose_touch,'k--')

np.savez( fileHandle, time=f['time'], trialFlag=f['trialFlag'],
        dFF=f['dFF'], dYY=f['dYY'], dRR=f['dRR'], ball=f['ball'], dlc=f['dlc'], 
        beh_labels=f['beh_labels'], stim=f['stim'], drink=f['drink'],
        dims=f['dims'], dims_in_um=f['dims_in_um'], im=f['im'], 
        scanRate=f['scanRate'], redTh=f['redTh'], grnTh=f['grnTh'],
        goodIds=f['goodIds'], oIsGood=f['oIsGood'], PIDdata=f['PIDdata'],
        aligned_centroids=f['aligned_centroids'], sucrose_touch=sucrose_touch ) 



