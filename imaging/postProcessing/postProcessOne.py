
import os
import sys
import matplotlib.pyplot as plt
from scipy import sparse, optimize
import numpy as np
import pickle
from scipy import io
# import smoothBehavior as beh
import source_process as sc

def postProcessOne(baseFolder, expNameHandle, savematfile):
    
    # obj = sc.scape(baseFolder)
    # obj.postProcess('F.mat', 'post_fromYcc', savematfile)
    obj = sc.scape(baseFolder)
    obj.postProcess(expNameHandle+'_raw.npz', expNameHandle, savematfile)
       

if __name__ == '__main__':
    mynargs = sys.argv
    if (len(mynargs)>1):
        print('processing folder:' + mynargs[1])
        rootFolder = mynargs[1]
        expID = mynargs[2]
    else:
        rootFolder = "/Users/evan/Dropbox/_sandbox/sourceExtraction/good/"
        expID = 'F_fromRed.mat'
    
    expNameHandle=expID.replace('/','_')
    savematfile=True
    postProcessOne(rootFolder, expNameHandle, savematfile)