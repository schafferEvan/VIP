#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Suite of functions to process scape data post-caiman

@author evan schaffer
"""
# Created Jan 3 2019


import sys
import numpy as np
from scipy import sparse, optimize, io
from skimage.restoration import denoise_tv_chambolle
import h5py


class smoothData:

    def __init__(self, baseFolder):
        self.baseFolder = baseFolder

    def loadMat7p3(self, file, structname, varname):
        f = h5py.File(self.baseFolder+'alignedBehavAndStim.mat')
        self.matVar = f[structname].__getitem__(structname).__getitem__(varname)
        self.matVar = np.transpose(self.matVar)


    def totVarSmoothData(self, data, weight):
        self.smoothData = np.zeros(np.shape(data))
        for i in range(0, np.shape(data)[0]):
            print(i)
            self.smoothData[i,:] = denoise_tv_chambolle(data[i,:], weight=weight)
        self.maxSm = np.amax(self.smoothData, axis=1)
        self.minSm = np.amin(self.smoothData, axis=1)

    def normalizeRawF(self, data):
        # explicitly passing in data rather than referencing attribute, so data can be 'Y' or 'C'
        self.scaled = np.zeros(np.shape(data))
        self.max = np.amax(data, axis=1)
        self.min = np.amin(data, axis=1)
        
        for i in range(0, np.shape(data)[0]):
            print(i)
            self.scaled[i,:] = (data[i,:]-self.min[i])/(self.max[i]-self.min[i])

    def getSmoothBeh(self):
        inputFile = 'alignedBehavAndStim.mat'
        outputFile = 'alignedBehavSmooth.mat'
        self.loadMat7p3(inputFile, '/alignedBehavior', 'legVar')
        print(np.shape(self.matVar))
        self.normalizeRawF(self.matVar)
        self.totVarSmoothData(self.scaled, 1.0)
        self.behSmooth = self.smoothData

        io.savemat(self.baseFolder+outputFile,{'behSmooth':self.behSmooth})


if __name__ == '__main__':
    mynargs = sys.argv
    if (len(mynargs)>1):
        print('processing folder:' + mynargs[1])
        baseFolder = mynargs[1]
    else:
        baseFolder = '/Users/evan/Dropbox/_sandbox/sourceExtraction/0824_f3r1/'
    
    obj = smoothData(baseFolder)
    obj.getSmoothBeh()
    
    
    
