
import os
import sys
import matplotlib.pyplot as plt
from scipy import sparse, optimize
import numpy as np
import pickle
from scipy import io
# import smoothBehavior as beh
import scape_caiman_postProcess as sc

def postProcessOne(baseFolder):
    
    obj = sc.scape(baseFolder)
    obj.postProcess('F.mat', 'post_fromYcc.mat')
    obj = sc.scape(baseFolder)
    obj.postProcess('F_fromRed.mat', 'post_fromRcc.mat')
       

if __name__ == '__main__':
    mynargs = sys.argv
    if (len(mynargs)>1):
        print('processing folder:' + mynargs[1])
        rootFolder = mynargs[1]
    else:
        rootFolder = "/Users/evan/Dropbox/_sandbox/sourceExtraction/good/"
    
    postProcessOne(rootFolder)