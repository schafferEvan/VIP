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
from scipy import sparse, optimize, io
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from skimage.restoration import denoise_tv_chambolle
from time import time
import scipy.cluster.hierarchy as sch
import colorsys
from sklearn.decomposition import FastICA
from numpy.polynomial.polynomial import Polynomial as poly
from scipy.ndimage.filters import gaussian_filter1d as gsm



class scape:

    def __init__(self, baseFolder):
        self.baseFolder = baseFolder
        self.A_tot = sparse.csc_matrix((1,1))
        self.C_tot = np.ndarray(shape=(1,1))
        self.YrA_tot = np.ndarray(shape=(1,1))
        self.b_tot = np.ndarray(shape=(1,1))
        self.f_tot = np.ndarray(shape=(1,1))
        self.Y = np.ndarray(shape=(1,1))
        self.Y0 = np.ndarray(shape=(1,1))
        self.Af = np.ndarray(shape=(1,1))
        self.b = np.ndarray(shape=(1,1))
        self.f = np.ndarray(shape=(1,1))
        self.tFit = np.ndarray(shape=(1,1))
        self.popt = np.ndarray(shape=(1,1))
        self.min = np.ndarray(shape=(1,1))
        self.max = np.ndarray(shape=(1,1))
        self.Yscaled = np.ndarray(shape=(1,1))
        self.fitCutPt = 0

    def loadMat(self, file, varname):
        self.mat=io.loadmat(self.baseFolder+file) 
        self.matVar=self.mat[varname]

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
        self.O = np.zeros(np.shape(self.YsmoothData))
        self.dOO = np.zeros(np.shape(self.YsmoothData))
        self.oIsGood = np.zeros((np.shape(self.YsmoothData)[0],1))
        nC = np.shape(self.O)[0]
        printFreq = int(nC/10)
        for i in range(0, nC):
            self.showProgress(i,printFreq)
            y = self.Ybl[i,:]
            r = self.Rbl[i,:]
            otmp = np.divide(y,r)
            otmp[np.flatnonzero(np.isinf(otmp))]=0
            otmp[np.flatnonzero(np.isnan(otmp))]=0
            self.O[i,:] = otmp

        self.Omax = np.amax(self.O, axis=1)
        self.Omin = np.amin(self.O, axis=1)

        for i in range(0, np.shape(self.O)[0]):
            dotmp = (self.O[i,:] - np.percentile(self.O[i,:], 10)) / np.percentile(self.O[i,:], 10)
            dotmp[np.flatnonzero(np.isinf(dotmp))]=0
            dotmp[np.flatnonzero(np.isnan(dotmp))]=0
            self.dOO[i,:] = dotmp
            if ((sum(self.O[i,:])>0) and (sum(self.dOO[i,:])>0)):
                self.oIsGood[i]=1
    

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



    def makePolyF0(self, data, ord):
        # explicitly passing in data rather than referencing attribute, so data can be 'Y' or 'C'
        self.tFit = np.linspace(0, np.shape(data)[1]-1, np.shape(data)[1])
        # self.popt = np.zeros((np.shape(data)[0],4))
        T = np.shape(data)[1]
        self.F0 = np.zeros(np.shape(data))
        self.rsq = np.zeros((np.shape(data)[0],1))

        for i in range(0, np.shape(data)[0]):
            print(i)
            f = poly.fit(self.tFit, data[i,:], ord)
            xx, self.F0[i,:] = f.linspace(T)
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

    def getGoodComponentsFull(self):
        # isGood = np.zeros((np.shape(self.Y)[0],1))
        ampTh = 500 #1500 #2000 # discard if max of trace is below this
        redTh = 100 #200 #2000 # discard if max of trace is below this
        magTh = 50 #2 #1  #discard if mean of dOO is greater than this (motion)
        minTh = 1 # discard if min is greater than this
        maxTh = 0.2 #0.3 # discard if max is smaller than this
        rgccTh = 0.98 #0.9 # discard units in which red and green are very correlated
        
        My = np.max(self.Y, axis=1)
        Mr = np.max(self.R, axis=1)
        Mo = np.max(self.dOO, axis=1)

        self.getDatacorr(self.dOO, self.Y)
        ogCorr = self.dataCorr
        self.getDatacorr(self.dOO, self.R)
        orCorr = self.dataCorr
        oMoreGreen = np.array(orCorr<ogCorr)
        self.oMoreGreen = oMoreGreen.flatten()

        self.ampIsGood = np.array(My>ampTh)
        self.redIsGood = np.array(Mr>redTh)
        rgccIsGood = np.array(self.rgCorr<rgccTh)
        self.rgccIsGood = rgccIsGood.flatten()
        self.minIsGood = np.array(np.min(self.dOO, axis=1)<minTh)
        self.maxIsGood = np.array(np.max(self.dOO, axis=1)>maxTh)
        self.magIsGood = np.array(np.mean(self.dOO, axis=1)<magTh)
        oIsGood = np.array(self.oIsGood>0)
        self.oIsGood = oIsGood.flatten()
        self.goodIds = self.ampIsGood & self.minIsGood & self.maxIsGood & self.magIsGood & self.rgccIsGood & self.oMoreGreen & self.redIsGood


    def getIdxList(self, longList, shortList):
        #self.trialFlag, self.trialFlagUnique
        idx = np.zeros(np.shape(shortList))
        for j in range(len(shortList)):
            idx[j] = np.where(longList==shortList[j])[0][0]
        return idx.astype(int)

    def showProgress(self, i, printFreq):
        if(not np.mod(i,printFreq)):
            print(i,end=" ")
            sys.stdout.flush()


    def saveSummary(self, filename):
        np.save(self.baseFolder+filename+'.npy',self.dOO)
        io.savemat(self.baseFolder+filename,{'Fsc':self.Yscaled,'Rsc':self.Rscaled,
            'F_max':self.Ymax,'F_min':self.Ymin,'R_max':self.Rmax,'R_min':self.Rmin,
            'Fsm':self.YsmoothData,'Rsm':self.RsmoothData,
            'dYY':self.dYY,'dRR':self.dRR,'dOO':self.dOO,
            'Ybl':self.Ybl,'Rbl':self.Rbl,
            'Ygoodsc':self.Ygoodsc,'Rgoodsc':self.Rgoodsc,
            'Y0sc':self.Y0sc,'R0sc':self.R0sc,
            'Fexp':self.Y0,'Rexp':self.R0,'O':self.O,
            'rsq':self.rsq,'oIsGood':self.oIsGood,'goodIds':self.goodIds,
            'ampIsGood':self.ampIsGood,'minIsGood':self.minIsGood,'maxIsGood':self.maxIsGood,'magIsGood':self.magIsGood,'rgccIsGood':self.rgccIsGood,'oMoreGreen':self.oMoreGreen,'redIsGood':self.redIsGood,
            'Ypopt':self.Ypopt,'Rpopt':self.Rpopt,
            })


    def importdata(self, inputFile):
        if inputFile.endswith('.mat'):
            self.loadMat(inputFile, 'FR')
            self.R = self.matVar
            self.loadMat(inputFile, 'F') 
            self.Y = self.matVar
            
            try:
                self.loadMat(inputFile, 'trialFlag')
                self.trialFlag = self.matVar
            except:
                # in old data with one trial, trialFlag is undefined
                self.trialFlag = np.ones(np.shape(self.Y)[1])
        elif inputFile.endswith('.npz'):
            d = np.load( inputFile )
            self.Y=d['Y']
            self.R=d['R']
            self.trialFlag = d['trialFlag']



  
    def postProcess(self, inputFile, outputFile):
        self.importdata(self.baseFolder+inputFile)

        print('\n normalizing red')
        self.normalizeRawF(self.R)
        self.Rscaled = self.scaled
        self.Rmax = self.max
        self.Rmin = self.min

        print('\n normalizing green')
        self.normalizeRawF(self.Y)
        self.Yscaled = self.scaled
        self.Ymax = self.max
        self.Ymin = self.min

        self.getDatacorr(self.R, self.Y)
        self.rgCorr = self.dataCorr

        try:
            self.loadMat(inputFile, 'trialFlag')
            self.trialFlag = self.matVar
        except:
            self.trialFlag = np.ones(np.shape(self.Y)[1])

        self.trialFlagMax = np.amax(self.trialFlag)

        # snrTh = -10**8 #15
        # expFitTh = 1 #0.99

        self.Rpopt = []
        self.Ypopt = []
        self.trialFlagUnique = np.unique(self.trialFlag)
        self.trList = self.getIdxList(self.trialFlag, self.trialFlagUnique)

        rdata = self.Rscaled
        ydata = self.Yscaled

        # do smoothing on data scaled from 0 to 1
        print('\n smoothing red data')
        self.totVarSmoothData(rdata, 1.0)
        self.RsmoothData = self.smoothData
        
        print('\n smoothing green data')
        self.totVarSmoothData(ydata, 1.0)
        self.YsmoothData = self.smoothData
        
        al = (0, 200)
        bnds = ((0, 1), (0.2, 2))
        # bnds = np.array( np.tile(al,(len(self.trList)+2,1)) )

        # self.makeQuantileDF0(RsmoothData, bnds, self.Rpopt)
        # self.makeQuantileF0_multiExp(RsmoothData, bnds, self.Rpopt)
        print('\n calculate red F0')
        self.makeQuantileDF0(self.RsmoothData, bnds, self.Rpopt)
        R0 = self.F0
        self.Rpopt = self.popt

        print('\n calculate green F0')
        self.makeQuantileDF0(self.YsmoothData, bnds, self.Ypopt)
        Y0 = self.F0
        self.Ypopt = self.popt

        # rescale smooth data to match scale of original data
        print('\n rescale red data')
        self.rescaleData(self.RsmoothData, self.Rmin, self.Rmax)
        Rgoodsc = self.rescaled

        print('\n rescale green data')
        self.rescaleData(self.YsmoothData, self.Ymin, self.Ymax)
        Ygoodsc = self.rescaled
        

        # rescale exponential fits to match scale of original data
        print('\n rescale red F0')
        self.rescaleData(R0, self.Rmin, self.Rmax)
        R0sc = self.rescaled

        print('\n rescale green F0')
        self.rescaleData(Y0, self.Ymin, self.Ymax)
        Y0sc = self.rescaled

        self.R0 = R0
        self.Y0 = Y0
        self.Rgoodsc = Rgoodsc
        self.Ygoodsc = Ygoodsc
        self.R0sc = R0sc
        self.Y0sc = Y0sc
        

        self.dYY = np.divide(self.Ygoodsc-self.Y0sc, self.Y0sc)
        self.dRR = np.divide(self.Rgoodsc-self.R0sc, self.R0sc)

        self.Ybl = np.transpose( np.transpose(self.Ygoodsc-self.Y0sc) + np.amin(self.Y0sc, axis=1) ) # bleach corrected but scaled like raw data
        self.Rbl = np.transpose( np.transpose(self.Rgoodsc-self.R0sc) + np.amin(self.R0sc, axis=1) )

        
        #self.dOO = np.divide(self.Ogoodsc-self.O0sc, self.O0sc)
        print('\n calculate O and dOO')
        self.make_O_and_dOO()
        self.getGoodComponentsFull()

        # dataToCluster = self.dOO[np.flatnonzero(self.goodIds),:]
        # self.computeCorr(dataToCluster)
        
        print('\n saving')
        self.saveSummary(outputFile)
        



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

    # obj = scape(baseFolder)
    # obj.postProcess('F.mat', 'post_fromYcc.mat')
    obj = scape(baseFolder)
    # obj.postProcess('F_fromRed.mat', 'post_fromRcc.mat')
    obj.postProcess('2019_06_26_Nsyb_NLS6s_walk_fly2_raw.npz', 'post_fromRcc.mat')



