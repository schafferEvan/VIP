#%%
import os
import sys
import glob
import csv
import numpy as np
import psutil
from scipy.ndimage.filters import gaussian_filter
import caiman as cm
from caiman.utils.visualization import nb_view_patches3d
import caiman.source_extraction.cnmf as cnmf
from caiman.components_evaluation import evaluate_components, estimate_components_quality_auto
from caiman.cluster import setup_cluster
from caiman.paths import caiman_datadir
import pickle
from scipy import io

from caiman.source_extraction.cnmf import cnmf as cnmf
from caiman.source_extraction.cnmf import params as params


#import matplotlib.pyplot as plt
#import bokeh.plotting as bpl
# stop the cluster if one exists
n_processes = psutil.cpu_count()
#print('using ' + str(n_processes) + ' processes')
#print("Stopping  cluster to avoid unnencessary use of memory....")
#sys.stdout.flush()  
#cm.stop_server()


baseFolder = '/mnt/ceph/users/evanschaffer/data/0110hunger_nobl/' #'/mnt/ceph/users/evanschaffer/data/'
mfilename = 'Yr_d1_110_d2_95_d3_249_order_C_frames_5701_.mmap' #'Yr_d1_124_d2_151_d3_277_order_C_frames_3791_.mmap' #'Yr_d1_112_d2_106_d3_257_order_C_frames_5064_.mmap' #'tYr_d1_80_d2_96_d3_241_order_C_frames_2823_.mmap' #'tiffYr_d1_144_d2_80_d3_274_order_C_frames_1005_.mmap' #'Yr_d1_131_d2_91_d3_255_order_C_frames_2703_.mmap' #'tiffsYr_d1_152_d2_126_d3_281_order_C_frames_871_.mmap' #'tiffsYr_d1_179_d2_156_d3_276_order_C_frames_3791_.mmap' #'tiffsYr_d1_114_d2_126_d3_281_order_C_frames_871_.mmap' #'reg20180817IRtestfly1Yr_d1_137_d2_73_d3_273_order_C_frames_8109_.mmap' #'reg20180801IRtestfly1Yr_d1_131_d2_91_d3_255_order_C_frames_8109_.mmap' #'reg20180816IRtestfly1Yr_d1_152_d2_73_d3_256_order_C_frames_8142_.mmap' #'reg20180801IRtestfly1Yr_d1_131_d2_91_d3_255_order_C_frames_8109_.mmap' #'reg20180817IRtestfly1Yr_d1_137_d2_73_d3_273_order_C_frames_8109_.mmap' #'reg20180816IRtestfly1Yr_d1_152_d2_73_d3_256_order_C_frames_8142_.mmap' #'reg20180817IRtestfly1Yr_d1_137_d2_73_d3_273_order_C_frames_8109_.mmap' #'reg20180816IRtestfly1Yr_d1_152_d2_73_d3_256_order_C_frames_8142_.mmap' #'sample301Yr_d1_131_d2_91_d3_255_order_C_frames_199_.mmap'
fname_new = [baseFolder + mfilename]

# Yr, dims, T = cm.load_memmap(fname_new)
# Y = np.reshape(Yr, dims + (T,), order='F')
# images = np.reshape(Yr.T, [T] + list(dims), order='F')    # reshape data in Python format (T x X x Y x Z)


# START CLUSTER
c, dview, n_processes = setup_cluster(
    backend='local', n_processes=24, single_thread=False)




# %% set up some parameters
# fnames = [os.path.join(caiman_datadir(), 'example_movies', 'demoMovie.tif')]
                        # file(s) to be analyzed
is_patches = True       # flag for processing in patches or not
fr = 8                  # approximate frame rate of data
decay_time = 50.0       # length of transient

if is_patches:          # PROCESS IN PATCHES AND THEN COMBINE
    rf = (10,10,40)             # half size of each patch
    stride = 10         # overlap between patches
    K = 30             # number of components in each patch
else:                   # PROCESS THE WHOLE FOV AT ONCE
    rf = None           # setting these parameters to None
    stride = None       # will run CNMF on the whole FOV
    K = 30              # number of neurons expected (in the whole FOV)

gSig = [3.0, 3.0, 3.0]  # expected half size of neurons
merge_thresh = 0.80     # merging threshold, max correlation allowed
p = 1                   # order of the autoregressive system
gnb = 1                 # global background order

params_dict = {'fnames': fname_new,
               'fr': fr,
               'decay_time': decay_time,
               'rf': rf,
               'stride': stride,
               'K': K,
               'gSig': gSig,
               'merge_thr': merge_thresh,
               'p': p,
               'nb': gnb,
               'low_rank_background': False}

opts = params.CNMFParams(params_dict=params_dict)
# %% Now RUN CaImAn Batch (CNMF)
cnm = cnmf.CNMF(n_processes, params=opts, dview=dview)
cnm = cnm.fit_file()
#cnm = cnm.fit(images)


# rf = (20,20,20) #(20,20,20) #(17,14,48) #(33, 28, 48)  # half-size of the patches in pixels. rf=25, patches are 50x50
# stride = (10,10,10) #(5,5,5) #(20, 10, 10)  # amounpl.it of overlap between the patches in pixels
# K = 200 #50 #250 #180  # number of neurons expected per patch
# gSig = [1.0, 2.0, 2.0] #[0.5, 1.0, 1.0]  # expected half size of neurons
# merge_thresh = 0.8  # merging threshold, max correlation allowed
# p = 1  # order of the autoregressive system
# save_results = False
# #%% RUN ALGORITHM ON PATCHES
# init_method = 'greedy_roi'
# alpha_snmf = None  # 10e2  # this controls sparsity
# images = np.reshape(Yr.T, [T] + list(dims), order='F')    # reshape data in Python format (T x X x Y x Z)
# cnm = cnmf.CNMF(n_processes, k=K, gSig=gSig, merge_thresh=0.8, p=p, dview=dview, Ain=None, rf=rf, stride=stride, memory_fact=2,
#                 method_init=init_method, alpha_snmf=alpha_snmf, only_init_patch=True, gnb=1, ssub=1, tsub=1, method_deconvolution='oasis')
# #skip_refinement=True

# cnm = cnm.fit(images)
# #cnm.save(baseFolder+'testOut')

A_tot = cnm.estimates.A
C_tot = cnm.estimates.C
YrA_tot = cnm.estimates.YrA
b_tot = cnm.estimates.b
f_tot = cnm.estimates.f
sn_tot = cnm.estimates.sn

print(('Number of components:' + str(A_tot.shape[-1])))

#with open(baseFolder+'C.csv', 'w+') as csvfile:
#    fw = csv.writer(csvfile, delimiter=',')
#    fw.writerows(C_tot)

#Af = A_tot.toarray() #note: this is needed because default format is scipy.sparse
#with open(baseFolder+'A.csv', 'w+') as csvfile:
#    fw = csv.writer(csvfile, delimiter=',')
#    fw.writerows(Af)

io.savemat(baseFolder+'all.mat',{'A':A_tot,'C':C_tot,'YrA':YrA_tot,'b':b_tot,'f':f_tot,'sn':sn_tot})

A_pkl = open(baseFolder+"A.pickle","wb")
pickle.dump(A_tot, A_pkl)
C_pkl = open(baseFolder+"C.pickle","wb")
pickle.dump(C_tot, C_pkl)
YrA_pkl = open(baseFolder+"YrA.pickle","wb")
pickle.dump(YrA_tot, YrA_pkl)
b_pkl = open(baseFolder+"b.pickle","wb")
pickle.dump(b_tot, b_pkl)
f_pkl = open(baseFolder+"f.pickle","wb")
pickle.dump(f_tot, f_pkl)
sn_pkl = open(baseFolder+"sn.pickle","wb")
pickle.dump(sn_tot, sn_pkl)



#Yr, dims, T = cm.load_memmap(fname_new[0])
#Y = np.reshape(Yr, dims + (T,), order='F')
#images = np.reshape(Yr.T, [T] + list(dims), order='F')    # reshape data in Python format (T x X x Y x Z)


#%% COMPONENT EVALUATION
# the components are evaluated in three ways:
#   a) the shape of each component must be correlated with the data
#   b) a minimum peak SNR is required over the length of a transient
#   c) each shape passes a CNN based classifier (this will pick up only neurons
#           and filter out active processes)

#min_SNR = 2.5       # peak SNR for accepted components (if above this, acept)
#rval_thr = 0.90     # space correlation threshold (if above this, accept)
#use_cnn = False     # use the CNN classifier
#min_cnn_thr = 0.95  # if cnn classifier predicts below this value, reject
#
#idx_components, idx_components_bad, SNR_comp, r_values, cnn_preds = \
#    estimate_components_quality_auto(images, cnm.estimates.A, cnm.estimates.C, cnm.estimates.b, cnm.estimates.f,
#                                     cnm.estimates.YrA, fr, decay_time, gSig, dims,
#                                     dview=dview, min_SNR=min_SNR,
#                                     r_values_min=rval_thr, use_cnn=use_cnn,
#                                     thresh_cnn_min=min_cnn_thr)

#idxcmp_pkl = open(baseFolder+"idx_components.pickle","wb")
#pickle.dump(idx_components, idxcmp_pkl)
#idx_components_bad_pkl = open(baseFolder+"idx_components_bad.pickle","wb")
#pickle.dump(idx_components_bad, idx_components_bad_pkl)
#SNR_comp_pkl = open(baseFolder+"SNR_comp.pickle","wb")
#pickle.dump(SNR_comp, SNR_comp_pkl)
#r_values_pkl = open(baseFolder+"r_values.pickle","wb")
#pickle.dump(r_values, r_values_pkl)

Amax = A_tot.max()
Amin = A_tot.min()

io.savemat(baseFolder+'all.mat',{'A':A_tot,'C':C_tot,'YrA':YrA_tot,'b':b_tot,'f':f_tot,'Amax':Amax,'Amin':Amin})

#%% STOP CLUSTER and clean up log files
#cm.stop_server()

#log_files = glob.glob('Yr*_LOG_*')
#for log_file in log_files:
#   os.remove(log_file)
