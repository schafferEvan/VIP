#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Suite of functions to process scape data after source extraction
@author evan schaffer
"""
# Created Jan 3 2019


import os
import sys
import glob
import numpy as np
import pickle
import warnings
from scipy import sparse, optimize, io
from scipy.optimize import minimize, nnls
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from skimage.restoration import denoise_tv_chambolle
from time import time
import scipy.cluster.hierarchy as sch
from sklearn.cluster import AgglomerativeClustering
import colorsys
from sklearn.decomposition import FastICA
from numpy.polynomial.polynomial import Polynomial as poly
from scipy.ndimage.filters import gaussian_filter1d as gsm
import pdb
import struct


class dataObj:
    def __init__(self):
        self.Y = np.ndarray(shape=(1,1))
        self.R = np.ndarray(shape=(1,1))
        self.trialFlag = np.ndarray(shape=(1,1))  
        self.time = np.ndarray(shape=(1,1))    
        self.ball = np.ndarray(shape=(1,1))  
        self.stim = np.ndarray(shape=(1,1))  
        self.drink = np.ndarray(shape=(1,1))     
        self.dlc = np.ndarray(shape=(1,1))  
        self.beh_labels = np.ndarray(shape=(1,1))   
        self.dims = np.ndarray(shape=(1,1))  
        self.dims_in_um = np.ndarray(shape=(1,1))  
        self.im = np.ndarray(shape=(1,1))   
        self.scanRate = np.ndarray(shape=(1,1))
        self.A = np.ndarray(shape=(1,1))      



class scape:

    def __init__(self, baseFolder):
        self.baseFolder = baseFolder
        self.raw = dataObj()
        self.good = dataObj()    

    def loadMat(self, file, varname):
        self.mat=io.loadmat(self.baseFolder+file) 
        self.matVar=self.mat[varname]

    def readPID(self, odorRun):
        self.PIDdata = {}
        # self.PIDdata['data'] = []
        # self.PIDdata['time'] = []
        self.PIDdata['odor_seq'] = []
        self.PIDdata['odor_initdelay'] = []
        self.PIDdata['odor_dur'] = []
        self.PIDdata['odor_isi'] = []
        self.PIDdata['odor_rep'] = []
        if os.path.isdir(self.baseFolder+'bin/'):
            infofile = glob.glob(self.baseFolder+'info/*run'+str(odorRun)+'*info.mat')[0]
            self.infodata=io.loadmat(infofile)
            self.PIDdata['odor_seq']=self.infodata['odor_seq']
            self.PIDdata['odor_initdelay'] = self.infodata['odor_inidelay']
            self.PIDdata['odor_dur'] = self.infodata['odor_duration']
            self.PIDdata['odor_isi'] = self.infodata['odor_isi']
            self.PIDdata['odor_rep'] = self.infodata['odor_rep']  

    def organizePIDRun(self):
        self.PIDdata['data'] = []
        self.PIDdata['time'] = []
        self.PIDdata['odorNum'] = []
        if self.raw.PIDAligned:
            self.PIDdata['odorNum'] = self.raw.PIDAligned['odorNum']
            self.PIDdata['data'] = self.raw.PIDAligned['data']
            self.PIDdata['time'] = self.raw.PIDAligned['time']

    def totVarSmoothData(self, data, weight):
        self.smoothData = np.zeros(np.shape(data))
        nC = np.shape(data)[0]
        printFreq = int(nC/10)
        for i in range(0, nC):
            self.showProgress(i,printFreq)
            self.smoothData[i,:] = denoise_tv_chambolle(data[i,:], weight=weight)
        self.maxSm = np.amax(self.smoothData, axis=1)
        self.minSm = np.amin(self.smoothData, axis=1)

    def normalizeRawF(self, data):
        # explicitly passing in data rather than referencing attribute, so data can be 'Y' or 'C'
        self.scaled = np.zeros(np.shape(data))
        self.max = np.amax(data, axis=1)
        self.min = np.amin(data, axis=1)
        
        nC = np.shape(data)[0]
        printFreq = int(nC/10)
        for i in range(0, nC):
            self.showProgress(i,printFreq)
            self.scaled[i,:] = (data[i,:]-self.min[i])/(self.max[i]-self.min[i])

    def rescaleData(self, data, origmin, origmax):
        # reinflate smooth data of chosen color to match scale of original data
        self.rescaled = np.zeros(np.shape(data))
        nC = np.shape(data)[0]
        printFreq = int(nC/10)
        for i in range(0, nC):
            self.showProgress(i,printFreq)
            sc = (origmax[i]-origmin[i]) #(origmax[i]-origmin[i])*(mymax[i]-mymin[i])
            ofs = origmin[i]             #origmin[i]+mymin[i]*(origmax[i]-origmin[i])
            self.rescaled[i,:] = data[i,:]*sc + ofs

    def getDatacorr(self, data1, data2):
        self.dataCorr = np.zeros((np.shape(data1)[0],1))
        
        for i in range(0, np.shape(data1)[0]):
            data = [data1[i,:],data2[i,:]]
            m = np.mean(data,1); m = m[:,np.newaxis]
            s = np.std(data,1); s = np.diag(np.reciprocal(s))
            X = np.matmul( s, (data - m) )
            cc = np.matmul(X, np.transpose(X))/np.shape(X)[1]
            self.dataCorr[i] = cc[0,1]


    def make_O_and_dOO(self):
        # compute ratio (O) and dO/O
        self.O = np.zeros(np.shape(self.Ygoodsc)) # YsmoothData
        self.dOO = np.zeros(np.shape(self.Ygoodsc))
        self.oIsGood = np.zeros((np.shape(self.Ygoodsc)[0],1))
        nC = np.shape(self.O)[0]
        printFreq = int(nC/10)
        for i in range(0, nC):
            self.showProgress(i,printFreq)
            y = self.Ygoodsc[i,:]
            r = self.Rgoodsc[i,:]
            otmp = np.divide(y,r)
            otmp[np.flatnonzero(np.isinf(otmp))]=0
            otmp[np.flatnonzero(np.isnan(otmp))]=0
            self.O[i,:] = otmp

        print('\n normalize O')
        self.normalizeRawF(self.O)
        self.Oscaled = self.scaled
        self.Omax = self.max
        self.Omin = self.min

        print('\n compute O0')
        popt = []
        bnds = ((0, 1), (0.2, 2))
        self.makeQuantileDF0(self.Oscaled, bnds, popt)
        self.O0 = self.F0
        self.Opopt = self.popt

        # rescale O0
        self.rescaleData(self.O0, self.Omin, self.Omax)
        self.O0sc = self.rescaled
        #compute dO/O
        self.dOO = np.divide(self.O-self.O0sc, self.O0sc)

        # pdb.set_trace()
        # below is legacy.  is any of it still needed?
        # self.Omax = np.amax(self.O, axis=1)
        # self.Omin = np.amin(self.O, axis=1)
        for i in range(self.O.shape[0]):
            # dotmp = (self.O[i,:] - np.percentile(self.O[i,:], 10)) / np.percentile(self.O[i,:], 10)
            # dotmp[np.flatnonzero(np.isinf(dotmp))]=0
            # dotmp[np.flatnonzero(np.isnan(dotmp))]=0
            # self.dOO[i,:] = dotmp
            if ((sum(self.O[i,:])>0) and (sum(self.dOO[i,:])>0)):
                self.oIsGood[i]=1
    

    def makeNNLS_F0_multiExp(self, data):
        # nonnegative lsq version 
        self.F0 = np.zeros(np.shape(data))
        self.Y0rescaled = np.zeros(np.shape(data))
        self.rsq = np.zeros((np.shape(data)[0],1))

        # make time/trail matrix
        ftime = np.arange(0,len(self.good.trialFlag),1)
        X = np.zeros((len(self.good.trialFlag),1+len(self.trList)))
        for i in range(len(self.trList)):
            if (i==len(self.trList)-1):
                X[self.trList[i]:,i] = 1
                X[self.trList[i]:,-1] = ftime[self.trList[i]:]-ftime[self.trList[i]]
            else:
                X[self.trList[i]:,i] = 1
                X[self.trList[i]:self.trList[i+1]-1,-1] = ftime[self.trList[i]:self.trList[i+1]-1]-ftime[self.trList[i]]
        X = -X

        nC = np.shape(data)[0]
        printFreq = int(nC/10)
        for i in range(0, nC):
            self.showProgress(i,printFreq)
            # pdb.set_trace()
            betaR,rnorm = nnls(X,np.log(data[i,:]))
            alphas = -np.cumsum(betaR[:-1])
            self.F0[i,:] = multiExpFun(ftime, np.exp(alphas), 0, self.trList, 1/betaR[-1])
            self.Y0rescaled[i,:] = np.multiply( self.F0[i,:], self.max[i]-self.min[i] ) + self.min[i]
            c = np.cov(self.F0[i,:], data[i,:])
            self.rsq[i] = (c[0,1]/np.sqrt(c[0,0]*c[1,1]))**2



    def makeQuantileF0_multiExp(self, data, bnds, poptPrev):
        # explicitly passing in data rather than referencing attribute, so data can be 'Y' or 'C'
        # self.tFit = np.linspace(0, np.shape(data)[1]-1, np.shape(data)[1])
        nobs = np.shape(data)[1]
        self.tFit = np.linspace(0, 1, nobs)
        time = np.arange(0,nobs,1)
        self.popt = np.zeros((np.shape(data)[0],len(bnds)))
        self.F0 = np.zeros(np.shape(data))
        self.Y0rescaled = np.zeros(np.shape(data))
        self.rsq = np.zeros((np.shape(data)[0],1))
        bnds = np.array(bnds)
        ics = np.ones((1,len(bnds)))

        for i in range(0, np.shape(data)[0]):
            print(i)
            res = minimize(multi_FobjFun, ics, [time,data[i,:],self.trList], method='SLSQP', bounds=bnds,
              options={'maxiter': 1000, 'ftol': 1e-8})
            self.popt[i,:] = res.x
            self.F0[i,:] = multiExpFun(time, res.x[:-2], res.x[-2], self.trList, res.x[-1])
            self.Y0rescaled[i,:] = np.multiply( self.F0[i,:], self.max[i]-self.min[i] ) + self.min[i]
            c = np.cov(self.F0[i,:], data[i,:])
            self.rsq[i] = (c[0,1]/np.sqrt(c[0,0]*c[1,1]))**2



    def makeQuantileDF0(self, data, bnds, poptPrev):
        # explicitly passing in data rather than referencing attribute, so data can be 'Y' or 'C'
        # self.tFit = np.linspace(0, np.shape(data)[1]-1, np.shape(data)[1])
        nobs = np.shape(data)[1]
        self.tFit = np.linspace(0, 1, nobs)
        self.dtFit = self.tFit[:-1]
        self.popt = np.zeros((np.shape(data)[0],3))
        self.F0 = np.zeros(np.shape(data))
        self.Y0rescaled = np.zeros(np.shape(data))
        self.rsq = np.zeros((np.shape(data)[0],1))
        # bnds = ((0, 1), (0, 2), (0.2, 2))
        bnds = np.array(bnds)

        nC = np.shape(data)[0]
        printFreq = int(nC/10)
        for i in range(0, nC):
            self.showProgress(i,printFreq)
            if (np.shape(poptPrev)[0]>0):
                bnds[1,:] = poptPrev[i,2] # if poptPrev is already populated, use existing value of tau
            ysm = gsm(data[i,:],200)
            dy = np.diff(ysm)*(nobs-1)
            res = minimize(dFobjFun, [0.5,1.0], [self.dtFit,dy], method='SLSQP', bounds=bnds)
            self.popt[i,0:3:2] = res.x
            x = np.argmin(data[i,:])
            self.popt[i,1] = data[i,x]-expFun(self.tFit[x], *self.popt[i,:])
            self.F0[i,:] = expFun(self.tFit, *self.popt[i,:])
            self.Y0rescaled[i,:] = np.multiply( self.F0[i,:], self.max[i]-self.min[i] ) + self.min[i]
            c = np.cov(self.F0[i,:], data[i,:])
            self.rsq[i] = (c[0,1]/np.sqrt(c[0,0]*c[1,1]))**2


    def makeQuantileF0(self, data, bnds):
        # explicitly passing in data rather than referencing attribute, so data can be 'Y' or 'C'
        # self.tFit = np.linspace(0, np.shape(data)[1]-1, np.shape(data)[1])
        self.tFit = np.linspace(0, 1, np.shape(data)[1])
        self.popt = np.zeros((np.shape(data)[0],3))
        self.F0 = np.zeros(np.shape(data))
        self.Y0rescaled = np.zeros(np.shape(data))
        self.rsq = np.zeros((np.shape(data)[0],1))
        # bnds = ((0, 1), (0, 2), (0.2, 2))

        for i in range(0, np.shape(data)[0]):
            print(i)
            res = minimize(FobjFun, [0.5,0.5,1], [self.tFit,data[i,:]], method='SLSQP', bounds=bnds)
            self.popt[i,:] = res.x
            self.F0[i,:] = expFun(self.tFit, *self.popt[i,:])
            self.Y0rescaled[i,:] = np.multiply( self.F0[i,:], self.max[i]-self.min[i] ) + self.min[i]
            c = np.cov(self.F0[i,:], data[i,:])
            self.rsq[i] = (c[0,1]/np.sqrt(c[0,0]*c[1,1]))**2



    def makeExpF0(self, data, bnds):
        # explicitly passing in data rather than referencing attribute, so data can be 'Y' or 'C'
        self.tFit = np.linspace(0, np.shape(data)[1]-1, np.shape(data)[1])
        self.popt = np.zeros((np.shape(data)[0],4))
        self.F0 = np.zeros(np.shape(data))
        self.Y0rescaled = np.zeros(np.shape(data))
        self.rsq = np.zeros((np.shape(data)[0],1))
        # self.dfF = np.zeros(np.shape(data))
        # bnds=( (0,0,-10), (10,0.2,10))
        
        for i in range(0, np.shape(data)[0]):
            print(i)
            
            self.popt[i,:], pcov = optimize.curve_fit(expFun, self.tFit, data[i,:], maxfev = 100000000, bounds=bnds )
            self.F0[i,:] = expFun(self.tFit, *self.popt[i,:])
            # self.popt[i,:], pcov = optimize.curve_fit(polyFun, self.tFit, data[i,:], maxfev = 100000000, bounds=bnds )
            # self.F0[i,:] = polyFun(self.tFit, *self.popt[i,:])
            self.Y0rescaled[i,:] = np.multiply( self.F0[i,:], self.max[i]-self.min[i] ) + self.min[i]
            #self.dfF[i,:] = np.divide( np.subtract( self.Yscaled[i,:], self.Y0[i,:] ), self.Y0[i,:] )
            #self.dfF[i,:] = self.dfF[i,:]-np.min(self.dfF[i,:])
            c = np.cov(self.F0[i,:], data[i,:])
            self.rsq[i] = (c[0,1]/np.sqrt(c[0,0]*c[1,1]))**2

    

    def fitExp(self, data, bnds):
        # explicitly passing in data rather than referencing attribute, so data can be 'Y' or 'C'
        self.tFit = np.linspace(0, np.shape(data)[1]-1, np.shape(data)[1])
        self.popt = np.zeros((np.shape(data)[0],4))
        self.rsq = np.zeros((np.shape(data)[0],1))
        self.max = np.amax(data, axis=1)
        self.min = np.amin(data, axis=1)

        #io.savemat('/Users/evan/Downloads/pix_popt.mat',{'popt':self.popt,'m':self.min,'M':self.max})
        for i in range(0, np.shape(data)[0]):
            print(i)
            try:
                dRe = (data[i,:]-self.min[i])/(self.max[i]-self.min[i])
                self.popt[i,:], pcov = optimize.curve_fit(expFun, self.tFit, dRe, maxfev = 100000000, bounds=bnds )
            except:
                print("failed to fit unit "+str(i))


    def computeCorr(self, data):
        print('computing correlation matrix')
        m = np.mean(data,1)
        m = m[:,np.newaxis]
        s = np.std(data,1)
        s = np.diag(np.reciprocal(s))
        X = np.matmul( s, (data - m) )
        self.cc = np.matmul(X, np.transpose(X))/np.shape(X)[1]


    def getGoodComponentsFull(self, redTh=None, grnTh=0):
        
        # grnTh = 25  # discard if max of trace is below this (not a cell). By design, this is a weak threshold to allow for inactive cells
        # redTh = 200 # discard if max of trace is below this (not a cell). This is the main threshold for accepting ROIs as cells
        # magTh = 50 #2 #1  #discard if mean of dOO is greater than this (motion)
        # minTh = 2 #1 # discard if min is greater than this (motion)
        # maxTh = 0.1 #0.2 # discard if max is smaller than this (just noise)
        # rgccTh = 0.95 #0.9 # discard units in which red and green are very correlated
        # motionTh = 10 # signal this large is probably motion artifact
        
        # My = np.max(self.good.Y, axis=1)
        # Mo = np.max(self.dOO, axis=1)

        # self.getDatacorr(self.dOO, self.good.Y)
        # ogCorr = self.dataCorr
        # self.getDatacorr(self.dOO, self.good.R)
        # orCorr = self.dataCorr
        # oMoreGreen = np.array(abs(orCorr)<abs(ogCorr))
        # self.oMoreGreen = oMoreGreen.flatten()

        if redTh is None:
            self.redIsGood = np.ones(self.good.Y.shape[0]) #np.array(Mr>redTh)
        else:
            Mr = np.max(self.good.R, axis=1)
            self.redIsGood = np.array(Mr>redTh)
        self.isNotMotion = np.ones(self.good.Y.shape[0]) #np.array(Mo<motionTh)
        self.ampIsGood = np.ones(self.good.Y.shape[0]) #np.array(My>grnTh)
        rgccIsGood = np.ones(self.good.Y.shape[0]) #np.array(self.rgCorr<rgccTh)
        self.rgccIsGood = rgccIsGood.flatten()
        # self.minIsGood = np.array(np.min(self.dOO, axis=1)<minTh)
        # self.maxIsGood = np.array(np.max(self.dOO, axis=1)>maxTh)
        # self.magIsGood = np.array(np.mean(self.dOO, axis=1)<magTh)
        oIsGood = np.array(self.oIsGood>0)
        self.oIsGood = oIsGood.flatten()
        self.goodIds = np.logical_and(self.oIsGood, self.redIsGood) #np.array(len(oIsGood)*[True]) #self.isNotMotion & self.ampIsGood & self.rgccIsGood & self.redIsGood
        # self.activeIds = self.maxIsGood
        #& self.minIsGood
        # pdb.set_trace()

        self.dOO = self.dOO[self.goodIds,:]
        self.dYY = self.dYY[self.goodIds,:]
        self.dRR = self.dRR[self.goodIds,:]
        self.good.A  = self.raw.A[:,self.goodIds] #self.good.A  = self.raw.A # 
        print('# Good Ids = '+str(self.goodIds.sum()))

    def hierCluster_fixedN(self,nClust):
        # version with prespecified cluster number
        cluster = AgglomerativeClustering(n_clusters=nClust, affinity='euclidean', linkage='ward')  
        cluster.fit_predict(self.dOO)  
        self.cluster_labels = cluster.labels_

    def hierCluster(self):
        # version using sklearn with variable cluster number
        # cluster = AgglomerativeClustering(
        #     n_clusters=None, affinity='cosine', 
        #     linkage='complete', distance_threshold=0.8)  #1.
        cluster = AgglomerativeClustering(
            n_clusters=None, affinity='euclidean', linkage='ward',distance_threshold=7)
        cluster.fit_predict(self.dOO)  
        # idx_new = np.argsort(cluster.labels_)
        self.cluster_labels = cluster.labels_
        print('found '+str(len(np.unique(cluster.labels_)))+' clusters')


    def getIdxList(self, longList, shortList):
        #self.trialFlag, self.trialFlagUnique
        idx = np.zeros(np.shape(shortList))
        for j in range(len(shortList)):
            idx[j] = np.where(longList==shortList[j])[0][0]
        return idx.astype(int)

    def showProgress(self, i, printFreq):
        if(not np.mod(i,printFreq)):
            print(i)
            sys.stdout.flush()

    def trimTrialStart(self,secsToTrim):
        # trim frames from beginning of each run
        extra_buffer = int(np.round( secsToTrim*self.raw.scanRate )) # number of frames to trim
        # self.good = self.raw

        isAkeeper = np.ones(np.shape( self.raw.trialFlag ))
        print('not trimming from first trial')
        for j in range(1,len( self.trialFlagUnique )):
            jbeg = np.where(self.raw.trialFlag ==self.trialFlagUnique[j])[0][0]
            isAkeeper[jbeg:(jbeg+extra_buffer)] = 0

        self.good.Y = self.raw.Y[:,isAkeeper[:,0]>0]
        self.good.R = self.raw.R[:,isAkeeper[:,0]>0]
        self.good.trialFlag = self.raw.trialFlag[isAkeeper[:,0]>0]
        # pdb.set_trace()
        self.good.time = self.raw.time[isAkeeper[:,0]>0]
        self.good.ball = self.raw.ball[isAkeeper[:,0]>0]
        if(self.raw.stim.shape[0]>2):
            self.good.stim = self.raw.stim[isAkeeper[:,0]>0]
            self.good.drink = self.raw.drink[isAkeeper[:,0]>0]
        else:
            self.good.stim = self.raw.stim
            self.good.drink = self.raw.drink
        if(self.raw.dlc.shape[0]>2):
            self.good.dlc = self.raw.dlc[:,isAkeeper[:,0]>0]
        else:
            self.good.dlc = self.raw.dlc
        self.good.dlc = self.good.dlc.T
        if(len(self.raw.beh_labels)>2):
            if self.raw.beh_labels.shape[1]>self.raw.beh_labels.shape[0]:
                self.raw.beh_labels = self.raw.beh_labels.T # in some datasets, saved raw labels were transposed
            
            self.good.beh_labels = self.raw.beh_labels[isAkeeper[:,0]>0]
        else:
            self.good.beh_labels = self.raw.beh_labels
        self.good.isAkeeper = isAkeeper



    def saveSummary(self, filename, savematfile):
        if savematfile:
            io.savemat(self.baseFolder+filename+'.mat',{'Fsc':self.Yscaled,'Rsc':self.Rscaled,
                'F_max':self.Ymax,'F_min':self.Ymin,'R_max':self.Rmax,'R_min':self.Rmin,
                'Fsm':self.YsmoothData,'Rsm':self.RsmoothData,
                'dYY':self.dYY,'dRR':self.dRR,'dOO':self.dOO,
                'Ygoodsc':self.Ygoodsc,'Rgoodsc':self.Rgoodsc,
                'Y0sc':self.Y0sc,'R0sc':self.R0sc,'O0sc':self.O0sc,
                'Fexp':self.Y0,'Rexp':self.R0,'O':self.O,
                'rsq':self.rsq,'oIsGood':self.oIsGood,'goodIds':self.goodIds,
                'ampIsGood':self.ampIsGood,'rgccIsGood':self.rgccIsGood, 'redTh':self.redTh, 'grnTh':self.grnTh,
                'redIsGood':self.redIsGood,'Ypopt':self.Ypopt,'Rpopt':self.Rpopt,
                }) # 'cluster_labels':self.cluster_labels,
            io.savemat(self.baseFolder+filename+'_Agood.mat',
                {'goodIds':self.goodIds, 'A':self.good.A, 'dims':self.raw.dims, 'centroids':self.raw.centroids, 'centroids_fromGreen':self.raw.centroids_fromGreen})

        np.savez( self.baseFolder+filename+'.npz', time=self.good.time, trialFlag=self.good.trialFlag,
                dFF=self.dOO, dYY=self.dYY, dRR=self.dRR, ball=self.good.ball, dlc=self.good.dlc, 
                beh_labels=self.good.beh_labels, stim=self.good.stim, drink=self.good.drink, 
                goodIds=self.goodIds, oIsGood=self.oIsGood, dims=self.raw.dims, dims_in_um=self.raw.dims_in_um, 
                im=self.raw.im, im_R=self.raw.im_R, scanRate=self.raw.scanRate, redTh=self.redTh, grnTh=self.grnTh, 
                aligned_centroids=[], PIDdata=self.PIDdata, isAkeeper=self.good.isAkeeper) 
        sparse.save_npz(self.baseFolder+filename+'_A.npz', self.good.A)


    def importdata(self, inputFile):
        if inputFile.endswith('.mat'):
            self.loadMat(inputFile, 'FR')
            self.raw.R = self.matVar
            self.loadMat(inputFile, 'F') 
            self.raw.Y = self.matVar
            
            try:
                self.loadMat(inputFile, 'trialFlag')
                self.raw.trialFlag = self.matVar
            except:
                # in old data with one trial, trialFlag is undefined
                self.raw.trialFlag = np.ones(np.shape(self.raw.Y)[1])
            self.raw.trialFlagMax = np.amax(self.raw.trialFlag)

        elif inputFile.endswith('.npz'):
            d = np.load( inputFile )
            self.raw.Y=d['Y']
            self.raw.R=d['R']
            self.raw.trialFlag = d['trialFlag']
            self.raw.time=d['time']
            self.raw.ball=d['ball']
            self.raw.dlc=d['dlc']
            self.raw.stim=d['stim']
            self.raw.drink=d['drink']
            self.raw.beh_labels=d['states'].T
            self.raw.dims=d['dims']
            self.raw.dims_in_um = d['tot_um_x'], d['tot_um_y'], d['tot_um_z']
            self.raw.im=d['im']
            if 'im_R' in d:
                self.raw.im_R=d['im_R']
            else:
                self.raw.im_R=None
            self.raw.scanRate=d['scanRate']
            self.raw.centroids=d['centroids']
            self.raw.centroids_fromGreen=d['centroids_fromGreen']
            self.raw.PIDAligned = d['PIDAligned']
            self.raw.A = sparse.load_npz( inputFile[:-7]+'A_raw.npz' )

        self.trialFlagUnique = np.unique(self.raw.trialFlag)
        self.trList = self.getIdxList(self.raw.trialFlag, self.trialFlagUnique)

        # fix normalization of R and Y (in extractF, these were calculated as sum over ROI instead of mean)
        Am = np.array(np.sum(self.raw.A,axis=0).T)
        # pdb.set_trace()
        self.raw.R = np.divide(self.raw.R, Am)
        self.raw.Y = np.divide(self.raw.Y, Am)





  
    def process(self, inputFile, outputFile, secsToTrim=10., savematfile=False, redTh=None, grnTh=0, odorRun=None):
        self.redTh = redTh
        self.grnTh = grnTh
        self.importdata(self.baseFolder+inputFile)
        self.trimTrialStart(secsToTrim)
        self.readPID(odorRun)

        print('\n normalizing red')
        self.normalizeRawF(self.good.R)
        self.Rscaled = self.scaled
        self.Rmax = self.max
        self.Rmin = self.min

        print('\n normalizing green')
        self.normalizeRawF(self.good.Y)
        self.Yscaled = self.scaled
        self.Ymax = self.max
        self.Ymin = self.min

        self.getDatacorr(self.good.R, self.good.Y)
        self.rgCorr = self.dataCorr

        # do smoothing on data scaled from 0 to 1
        print('\n smoothing red data')
        self.totVarSmoothData(self.Rscaled, 1.0)
        self.RsmoothData = self.smoothData
        
        print('\n smoothing green data')
        self.totVarSmoothData(self.Yscaled, 1.0)
        self.YsmoothData = self.smoothData
        
        self.Rpopt = []
        self.Ypopt = []
        al = (0, 200)
        bnds = ((0, 1), (0.2, 2))
        # bnds = np.array( np.tile(al,(len(self.trList)+2,1)) )

        # self.makeQuantileDF0(RsmoothData, bnds, self.Rpopt)
        # self.makeQuantileF0_multiExp(RsmoothData, bnds, self.Rpopt)
        print('\n calculate red F0')
        # self.makeNNLS_F0_multiExp(self.RsmoothData)
        self.makeQuantileDF0(self.RsmoothData, bnds, self.Rpopt)
        self.R0 = self.F0
        self.Rpopt = self.popt

        print('\n calculate green F0')
        # self.makeNNLS_F0_multiExp(self.YsmoothData)
        self.makeQuantileDF0(self.YsmoothData, bnds, self.Ypopt)
        self.Y0 = self.F0
        self.Ypopt = self.popt

        # rescale smooth data to match scale of original data
        print('\n rescale red data')
        self.rescaleData(self.RsmoothData, self.Rmin, self.Rmax)
        self.Rgoodsc = self.rescaled

        print('\n rescale green data')
        self.rescaleData(self.YsmoothData, self.Ymin, self.Ymax)
        self.Ygoodsc = self.rescaled
        

        # rescale exponential fits to match scale of original data
        print('\n rescale red F0')
        self.rescaleData(self.R0, self.Rmin, self.Rmax)
        self.R0sc = self.rescaled

        print('\n rescale green F0')
        self.rescaleData(self.Y0, self.Ymin, self.Ymax)
        self.Y0sc = self.rescaled

        self.dYY = np.divide(self.Ygoodsc-self.Y0sc, self.Y0sc)
        self.dRR = np.divide(self.Rgoodsc-self.R0sc, self.R0sc)

        # self.Ybl = np.transpose( np.transpose(self.Ygoodsc-self.Y0sc) + np.amin(self.Y0sc, axis=1) ) # bleach corrected but scaled like raw data
        # self.Rbl = np.transpose( np.transpose(self.Rgoodsc-self.R0sc) + np.amin(self.R0sc, axis=1) )

        
        #self.dOO = np.divide(self.Ogoodsc-self.O0sc, self.O0sc)
        print('\n calculate O and dOO (ratiometric dFF)')
        self.make_O_and_dOO()

        print('\n find and remove bad units')
        self.getGoodComponentsFull(redTh, grnTh)

        # dataToCluster = self.dOO[np.flatnonzero(self.goodIds),:]
        # self.computeCorr(dataToCluster)
        # if savematfile:
        #     print('clustering')
        #     self.hierCluster() #_fixedN(20)
        
        print('\n saving')
        self.saveSummary(outputFile, savematfile)
        



def expFun(t, a, b, tau):
    return b + a*np.exp(- t/tau)

def dExpFun(t, a, tau):
    return  - a/tau*np.exp(- t/tau)

def multiExpFun(time, a, b, idx, tau):
    x = b + np.zeros(len(time))
    for i in range(len(idx)):
        T = idx[i]
        x[T:] = a[i]*np.exp(- (time[T:]-T)/tau)
    return x

def dFobjFun(params, data):
    # this is the objective function to optimize in quantile regression, with 100*q the quantile in %
    q = 0.5
    a, tau = params
    [X,Y] = data
    Yhat = dExpFun( X,a,tau)
    h = (Yhat>Y)
    hn = np.flatnonzero(np.logical_not(h))
    hp = np.flatnonzero(h)
    return (q-1)*np.sum( Y[hp] - Yhat[hp] ) + q*np.sum( Y[hn] - Yhat[hn] )

def FobjFun(params, data):
    # this is the objective function to optimize in quantile regression, with 100*q the quantile in %
    q = 0.0001
    a, b, tau = params
    [X,Y] = data
    Yhat = expFun( X,a,b,tau)
    h = (Yhat>Y)
    hn = np.flatnonzero(np.logical_not(h))
    hp = np.flatnonzero(h)
    return (q-1)*np.sum( Y[hp] - Yhat[hp] ) + q*np.sum( Y[hn] - Yhat[hn] )

def multi_FobjFun(params, data):
    # this is the objective function to optimize in quantile regression, with 100*q the quantile in %
    q = 0.5
    a = params[:-2]
    b = params[-2]
    tau = params[-1]
    [X,Y,trList] = data
    
    Yhat = multiExpFun(X, a, b, trList, tau)
    h = (Yhat>Y)
    hn = np.flatnonzero(np.logical_not(h))
    hp = np.flatnonzero(h)
    objOut = (q-1)*np.sum( Y[hp] - Yhat[hp] ) + q*np.sum( Y[hn] - Yhat[hn] )
    return objOut



if __name__ == '__main__':

    baseFolder = '/Volumes/SCAPEdata1/finalData/2019_06_26_Nsyb_NLS6s_walk/fly2/Yproj/' #'/Volumes/dataFast/sample/2019_06_26_Nsyb_NLS6s_walk/fly2/Yproj/'

    savematfile = False #True
    secsToTrim=10.
    # obj = scape(baseFolder)
    # obj.postProcess('F.mat', 'post_fromYcc.mat')
    obj = scape(baseFolder)
    # obj.postProcess('F_fromRed.mat', 'post_fromRcc.mat')
    obj.process('2019_06_26_Nsyb_NLS6s_walk_fly2_raw.npz', '2019_06_26_Nsyb_NLS6s_walk_fly2.npz',secsToTrim, savematfile)




