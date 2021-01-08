# Scripts for running post-processing on Terremoto


import os
import sys
import matplotlib.pyplot as plt
from scipy import sparse, optimize
import numpy as np
import pickle
from scipy import io
# import smoothBehavior as beh
import source_process as sc


def postProcessOne(baseFolder, expNameHandle, secsToTrim, savematfile, redTh, grnTh, odorNum):
    
    # obj = sc.scape(baseFolder)
    # obj.postProcess('F.mat', 'post_fromYcc', savematfile)
    obj = sc.scape(baseFolder)
    obj.process(expNameHandle+'_raw.npz', expNameHandle, secsToTrim, savematfile, redTh, grnTh, odorNum)
       

if __name__ == '__main__':
    mynargs = sys.argv
    if (len(mynargs)>1):
        print('processing folder:' + mynargs[1])
        rootFolder = mynargs[1]
        expID = mynargs[2]
        if (len(mynargs)>3):
            redTh=mynargs[3]
            redTh=float(redTh)
        else:
            redTh=100
        if (len(mynargs)>4):
            grnTh=float(mynargs[4])
        else:
            grnTh=0
        if (len(mynargs)>5):
            secsToTrim=float(mynargs[5])
        else:
            secsToTrim=20.
        if (len(mynargs)>6):
            odorNum=float(mynargs[6])
        else:
            odorNum=None
    else:
        expID = '2018_08_24_fly2/fly2'
        rootFolder = '/moto/axs/projects/registered/'+expID+'/Yproj/'
        secsToTrim = 0
        redTh = 100 
        grnTh = 0
        odorNum = None  
    
    expNameHandle=expID.replace('/','_')
    savematfile=True
    postProcessOne(rootFolder, expNameHandle, secsToTrim, savematfile, redTh, grnTh, odorNum)