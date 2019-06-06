#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Suite of functions to process scape data post-caiman

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




# Yr, dims, T = cm.load_memmap(fname_new)
# Y = np.reshape(Yr, dims + (T,), order='F')
# images = np.reshape(Yr.T, [T] + list(dims), order='F')    # reshape data in Python format (T x X x Y x Z)

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

    def makeRawF(self):
        print('making raw F')
        self.Y = np.zeros(np.shape(self.C_tot))
        self.Af = self.A_tot #.toarray()
        print('1')
        
        self.Amsk = np.zeros(np.shape(self.A_tot), dtype=bool)    
        self.Amsk[self.A_tot.nonzero()]=True #self.Af>0

        print('2')
        self.Asum = np.sum(self.Af,0)
        print('3')
        self.bmsk = np.zeros(np.shape(self.Asum))
        print('starting loop')
        print(np.shape(self.bmsk))
        print(np.shape(self.Amsk))
        print(np.shape(self.b))
        print(np.shape(self.Amsk[:,0]))
        t = time()
        
        self.bmsk = np.matmul( np.asarray(self.b), self.Amsk)
        # for i in range(0, np.shape(self.C_tot)[0]):
        #     self.bmsk[0,i] = np.sum(self.b[self.Amsk[:,i]])
        
        elapsed = time() - t; print(elapsed)
        self.Y = np.matmul(self.Asum, self.C_tot+self.YrA_tot ) + np.outer(self.bmsk,self.f)
        elapsed = time() - t; print(elapsed)

    def makeRawF_loop(self):
        print('making raw F')
        self.Y = np.zeros(np.shape(self.C_tot))
        self.Af = self.A_tot.toarray()
        
        for i in range(0, np.shape(self.C_tot)[0]):
            print(i)
            self.Amsk = self.Af[:,i]>0
            self.Asum = np.sum(self.Af[:,i])
            self.bmsk = np.sum(self.b[self.Amsk])
            self.Y[i,:] = self.Asum*(self.C_tot[i,:]+self.YrA_tot[i,:]) + self.bmsk*self.f
        #scape_pkl = open(self.baseFolder+"Y.pickle","wb")
        #pickle.dump(self.Y, scape_pkl)

    def loadRawF(self):
        self.Y = pickle.load( open(self.baseFolder+"Y.pickle","rb") )

    def loadCaiman(self):
        self.A_tot = pickle.load( open(self.baseFolder+"A.pickle","rb") )
        self.C_tot = pickle.load( open(self.baseFolder+"C.pickle","rb") )
        self.YrA_tot = pickle.load( open(self.baseFolder+"YrA.pickle","rb") )
        self.b_tot = pickle.load( open(self.baseFolder+"b.pickle","rb") )
        self.f_tot = pickle.load( open(self.baseFolder+"f.pickle","rb") )
        self.b = np.sum(self.b_tot,axis=1)
        self.f = np.sum(self.f_tot,axis=0)

    def loadMat(self, file, varname):
        self.mat=io.loadmat(self.baseFolder+file) 
        self.matVar=self.mat[varname]

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

    def rescaleData(self, data, origmin, origmax):
        # reinflate smooth data of chosen color to match scale of original data
        self.rescaled = np.zeros(np.shape(data))
        
        for i in range(0, np.shape(data)[0]):
            print(i)
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

        for i in range(0, np.shape(self.O)[0]):
            print(i)
            # y = self.YsmoothData[i,:]*(self.Ymax[i]-self.Ymin[i])+self.Ymin[i]
            # r = self.RsmoothData[i,:]*(self.Rmax[i]-self.Rmin[i])+self.Rmin[i]
            y = self.Ybl[i,:]
            r = self.Rbl[i,:]
            otmp = np.divide(y,r)
            otmp[np.flatnonzero(np.isinf(otmp))]=0
            otmp[np.flatnonzero(np.isnan(otmp))]=0
            self.O[i,:] = otmp
            # m = min(self.O[i,:])
            # dotmp = (self.O[i,:]-m)/m
            # dotmp[np.flatnonzero(np.isinf(dotmp))]=0
            # dotmp[np.flatnonzero(np.isnan(dotmp))]=0
            # self.dOO[i,:] = dotmp

        # bnds = ((0, 1), (0.2, 2))
        # self.makeQuantileDF0(self.O, bnds)
        # self.O0 = self.F0
        self.Omax = np.amax(self.O, axis=1)
        self.Omin = np.amin(self.O, axis=1)
        # self.rescaleData(self.O0, self.Omin, self.Omax)
        # self.O0sc = self.rescaled

        for i in range(0, np.shape(self.O)[0]):
            # dotmp = np.divide( np.subtract( self.O[i,:], self.O0[i,:] ), self.O0[i,:] )
            # m = min(dotmp)
            # dotmp = dotmp-m
            dotmp = (self.O[i,:] - np.percentile(self.O[i,:], 10)) / np.percentile(self.O[i,:], 10)
            dotmp[np.flatnonzero(np.isinf(dotmp))]=0
            dotmp[np.flatnonzero(np.isnan(dotmp))]=0
            self.dOO[i,:] = dotmp
            if ((sum(self.O[i,:])>0) and (sum(self.dOO[i,:])>0)):
                self.oIsGood[i]=1
    

    def calculate_I_and_dII(self):
        # Leifer lab method of separating signal (I) from motion artifacts using ICA
        self.I = np.zeros(np.shape(self.YsmoothData))
        self.dII = np.zeros(np.shape(self.YsmoothData))
        self.icflipped = np.zeros((np.shape(self.YsmoothData)[0],1))
        self.sigIc = np.zeros((np.shape(self.YsmoothData)[0],1))
        ica = FastICA(n_components=2)

        print(' PROBLEM HERE IS THAT ICA OUTPUT IS WHITENED, SO I LOOKS OK, BUT DII IS WRONG BECAUSE OF DENOMINATOR')
        for i in range(0, np.shape(self.I)[0]):
            #if (self.oIsGood[i]>0):
            try:
                data = [self.YsmoothData[i,:],self.RsmoothData[i,:]]
                m = np.mean(data,1); m = m[:,np.newaxis]
                s = np.std(data,1); s = np.diag(np.reciprocal(s))
                X = np.matmul( s, (data - m) )
                ic = ica.fit_transform(np.transpose(data))  # Reconstruct signals
                # ic = ica.fit_transform(np.transpose(X),whiten=False)  # Reconstruct signals
                s_ = np.transpose(ic)
                mm = np.mean(s_,1); mm = mm[:,np.newaxis]
                ss = np.std(s_,1); ss = np.diag(np.reciprocal(ss))
                Xic = np.matmul( ss, (s_ - mm) )
                cc = np.matmul(Xic, np.transpose(X))/np.shape(X)[1] #
                # sigcom = np.argmin(np.abs(cc[:,1])) # signal component is the one less correlated with red channel (ignoring sign)
                sigcom = np.argmax(np.abs(cc[:,0])) # signal component is the one more correlated with green channel (ignoring sign)
                self.sigIc[i] = cc[sigcom,0]
                if (cc[sigcom,0]<0): 
                    self.I[i,:] = -s_[sigcom,:] #np.matmul(-s_[sigcom,:],1/s[sigcom,sigcom]) + m[sigcom]
                    self.icflipped[i] = 1
                else:
                    self.I[i,:] = s_[sigcom,:] #np.matmul(s_[sigcom,:],1/s[sigcom,sigcom]) + m[sigcom]
            except:
                print('unable to perform ICA on unit')

        # bnds = ((0, 1), (0.2, 2))
        # self.makeQuantileDF0(self.I, bnds)
        # self.I0 = self.F0
        # for i in range(0, np.shape(self.I)[0]):
        #     try:
        #         self.dII[i,:] = np.divide( np.subtract( self.I[i,:], self.I0[i,:] ), self.I0[i,:] )
        #     except:
        #         print('unable to calculate dII on unit')
            

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

        for i in range(0, np.shape(data)[0]):
            print(i)
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

            # try:
            #     self.popt[i,:], pcov = optimize.curve_fit(expFun, self.tFit, data[i,:], maxfev = 100000000, bounds=bnds )
            #     self.F0[i,:] = expFun(self.tFit, *self.popt[i,:])
            #     self.Y0rescaled[i,:] = np.multiply( self.F0[i,:], self.max[i]-self.min[i] ) + self.min[i]
            #     #self.dfF[i,:] = np.divide( np.subtract( self.Yscaled[i,:], self.Y0[i,:] ), self.Y0[i,:] )
            #     #self.dfF[i,:] = self.dfF[i,:]-np.min(self.dfF[i,:])
            #     c = np.cov(self.F0[i,:], data[i,:])
            #     self.rsq[i] = (c[0,1]/np.sqrt(c[0,0]*c[1,1]))**2
            # except:
            #     print("failed to fit unit "+str(i))

            #bounds=( (0,0,-1e4), (1e5,1e4,1e7))
            #self.Y0[i,:] = 
        #scape_pkl = open(self.baseFolder+"popt.pickle","wb")
        #pickle.dump(self.popt, scape_pkl)

    

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
        #io.savemat('/Users/evan/Downloads/pix_popt.mat',{'popt':self.popt,'m':self.min,'M':self.max})


    def getSNR(self, rawData, smoothData):
        self.snr = np.zeros((np.shape(rawData)[0],1))
        for i in range(0, np.shape(rawData)[0]):
            self.snr[i] = np.var(smoothData[i,:])/np.var(rawData[i,:]-smoothData[i,:])

    def computeCorr(self, data):
        print('computing correlation matrix')
        m = np.mean(data,1)
        m = m[:,np.newaxis]
        s = np.std(data,1)
        s = np.diag(np.reciprocal(s))
        X = np.matmul( s, (data - m) )
        self.cc = np.matmul(X, np.transpose(X))/np.shape(X)[1]

    def hierCluster(self, corrData):
        print('computing hierarchical clustering')
        d = sch.distance.pdist(corrData)   # vector of ('55' choose 2) pairwise distances
        L = sch.linkage(d, method='weighted') #'complete' #average
        # self.clustInd = sch.fcluster(L, 0.1*d.max(), 'distance')
        # self.clustInd = sch.fcluster(L, 0.5, 'inconsistent')
        self.clustInd = sch.fcluster(L, 3.0, 'distance') #4.0 #0.5*d.max() # 4 too low.  6 too high.

    def plotClusters(self, data):
        # plt.figure(0)
        s = list(set(self.clustInd))
        self.clustData = np.zeros((len(s),np.shape(data)[1]))
        npl = np.ceil(np.sqrt(len(s)))
        time = np.linspace(0, np.shape(data)[1]-1, np.shape(data)[1])
        for j in range(0, len(s)):
            k = [i for i,x in enumerate(s) if x==s[j]] # get elements equal to s[j]
            # ax=plt.subplot(npl,npl,j+1)
            idata = data[np.flatnonzero(self.clustInd==(j+1)),:]
            M = np.max(idata,0)
            m = np.min(idata,0)
            # plt.plot(m,color=(.5,.5,.5))
            # plt.plot(M,color=(.5,.5,.5))
            mn = np.mean(idata,0)
            self.clustData[j,:] = mn
        #     plt.plot(mn, color=colorsys.hsv_to_rgb(j/len(s),1,1))
        #     plt.axis('off')
        # plt.show()


    def computeICA(self, data, nClust):
        # Compute ICA
        ica = FastICA(n_components=nClust)
        self.S_ = ica.fit_transform(np.transpose(data))  # Reconstruct signals
        self.A_ = ica.mixing_  # Get estimated mixing matrix
        makeICAfig=False
        if makeICAfig:
            plt.figure()
            cmap = get_cmap('tab10') #Dark2
            colors =  cmap(np.linspace(0, nClust-1, nClust)/nClust)
            plt.title('ICA recovered signals')
            for j in range(0,nClust):
                plt.plot(self.S_[:,j], color=colors[j])
            plt.show()


    def getGoodComponents(self, snrTh, expFitTh, data):
        self.isGood = (self.snr>snrTh)
        self.goodData = data[np.flatnonzero(self.isGood),:]

    def getGoodComponentsFull(self):
        isGood = np.zeros((np.shape(self.Y)[0],1))
        ampTh = 1500 #2000 # discard if max of trace is below this
        redTh = 200 #2000 # discard if max of trace is below this
        magTh = 2 #1  #discard if mean of dOO is greater than this (motion)
        minTh = 1 # discard if min is greater than this
        maxTh = 0.2 #0.3 # discard if max is smaller than this
        rgccTh = 0.9 # discard units in which red and green are very correlated
        
        My = np.max(self.Y, axis=1)
        Mr = np.max(self.R, axis=1)
        Mo = np.max(self.dOO, axis=1)

        self.getDatacorr(self.dOO, self.Y)
        ogCorr = self.dataCorr
        self.getDatacorr(self.dOO, self.R)
        orCorr = self.dataCorr
        oMoreGreen = np.array(orCorr<ogCorr)
        oMoreGreen = oMoreGreen.flatten()

        ampIsGood = np.array(My>ampTh)
        redIsGood = np.array(Mr>redTh)
        rgccIsGood = np.array(self.rgCorr<rgccTh)
        rgccIsGood = rgccIsGood.flatten()
        minIsGood = np.array(np.min(self.dOO, axis=1)<minTh)
        maxIsGood = np.array(np.max(self.dOO, axis=1)>maxTh)
        magIsGood = np.array(np.mean(self.dOO, axis=1)<magTh)
        oIsGood = np.array(self.oIsGood>0)
        oIsGood = oIsGood.flatten()
        self.goodIds = ampIsGood & oIsGood & minIsGood & maxIsGood & magIsGood & rgccIsGood & oMoreGreen & redIsGood



    def showF0hists(self):
        plt.figure(1)
        plt.subplot(311)
        plt.hist(self.popt[:,0], 200, density=True, facecolor='r', alpha=0.75)
        plt.subplot(312)
        plt.hist(self.popt[:,1], 200, density=True, facecolor='g', alpha=0.75)
        plt.subplot(313)
        plt.hist(self.popt[:,2], 200, density=True, facecolor='b', alpha=0.75)
        plt.show()

    def showFitSample(self, idx):
        plt.plot(self.tFit, self.YsmoothData[idx,:], 'b-', label='data')
        # plt.plot(self.tFit, expFun(self.tFit, *self.popt[idx,:]), 'r--', label='fit')
        # plt.plot(self.tFit, polyFun(self.tFit, *self.popt[idx,:]), 'r--', label='fit')
        plt.plot(self.tFit, self.Y0[idx,:], 'r--', label='fit')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.legend()
        plt.show()

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
            'snr':self.snr,'isGood':self.isGood,'rsq':self.rsq,'oIsGood':self.oIsGood,'goodIds':self.goodIds,
            'Ypopt':self.Ypopt,'Rpopt':self.Rpopt,'rgCorr':self.rgCorr,
            'I':self.I,'icflipped':self.icflipped,'sigIc':self.sigIc,
            'cc':self.cc,'clustInd':self.clustInd,'clustData':self.clustData,'ICA':self.S_})

    #'Osc':self.Oscaled,'Osm':self.OsmoothData,'OgoodData':self.OgoodData,'Ogoodsc':self.Ogoodsc,'O0sc':self.O0sc,'Oexp':self.O0,
    #'F_maxSm':self.YmaxSm,'F_minSm':self.YminSm,'R_maxSm':self.RmaxSm,'R_minSm':self.RminSm,
    #'O0sc':self.O0sc,'Oexp':self.O0,'dII':self.dII,'Iexp':self.I0,

    def postProcess(self, inputFile, outputFile):
        self.loadMat(inputFile, 'FR')
        self.R = self.matVar
        self.normalizeRawF(self.R)
        self.Rscaled = self.scaled
        self.Rmax = self.max
        self.Rmin = self.min

        self.loadMat(inputFile, 'F') # do Y after R so other attributes aren't overwritten
        self.Y = self.matVar
        self.normalizeRawF(self.Y)
        self.Yscaled = self.scaled
        self.Ymax = self.max
        self.Ymin = self.min

        self.getDatacorr(self.R, self.Y)
        self.rgCorr = self.dataCorr
        # self.O = np.divide(self.Y,self.R)
        # self.normalizeRawF(self.O)
        # self.Oscaled = self.scaled
        # self.Omax = self.max
        # self.Omin = self.min
        # ****************************************
        # ****************************************

        
        # (V) -try loadMat -> trialFLag, catch trialFLag=ones
        try:
            self.loadMat(inputFile, 'trialFlag')
            self.trialFlag = self.matVar
        except:
            self.trialFlag = np.ones(np.shape(self.Y)[1])

        # (V) -get max of trialFlag
        self.trialFlagMax = np.amax(self.trialFlag)
        # plt.plot(self.trialFlag)
        # plt.show()
        # -put most of what follows in loop over elements in trialFlag to process individually and concatenate
        # -for exp fit, after first cycle through loop (trial>1), keep tau fixed and just optimize other parameters
        
        snrTh = -10**8 #15
        expFitTh = 1 #0.99
        bnds = ((0, 1), (0.2, 2))
        self.Rpopt = []
        self.Ypopt = []
        print('**********')
        print(np.int(1+np.round(self.trialFlagMax)))
        print('**********')
        for i in range(1, np.int(1+np.round(self.trialFlagMax))):
            rdata = self.Rscaled[:,np.flatnonzero(self.trialFlag==(i))]
            ydata = self.Yscaled[:,np.flatnonzero(self.trialFlag==(i))]

            # do smoothing on data scaled from 0 to 1
            self.totVarSmoothData(rdata, 1.0)
            RsmoothData = self.smoothData
            # self.RmaxSm = self.maxSm
            # self.RminSm = self.minSm
            self.totVarSmoothData(ydata, 1.0)
            YsmoothData = self.smoothData
            # self.YmaxSm = self.maxSm
            # self.YminSm = self.minSm
            
            # compute exponential fit on smoothed data
            self.makeQuantileDF0(RsmoothData, bnds, self.Rpopt)
            R0 = self.F0
            self.Rpopt = self.popt
            self.makeQuantileDF0(YsmoothData, bnds, self.Ypopt)
            Y0 = self.F0
            self.Ypopt = self.popt

            # rescale smooth data to match scale of original data
            self.rescaleData(YsmoothData, self.Ymin, self.Ymax)
            Ygoodsc = self.rescaled
            self.rescaleData(RsmoothData, self.Rmin, self.Rmax)
            Rgoodsc = self.rescaled

            # rescale exponential fits to match scale of original data
            self.rescaleData(Y0, self.Ymin, self.Ymax)
            Y0sc = self.rescaled
            self.rescaleData(R0, self.Rmin, self.Rmax)
            R0sc = self.rescaled

            if (i==1):
                self.RsmoothData = RsmoothData
                self.YsmoothData = YsmoothData
                self.R0 = R0
                self.Y0 = Y0
                self.Rgoodsc = Rgoodsc
                self.Ygoodsc = Ygoodsc
                self.R0sc = R0sc
                self.Y0sc = Y0sc
            else:
                self.RsmoothData = np.concatenate((self.RsmoothData, RsmoothData), axis=1)
                self.YsmoothData = np.concatenate((self.YsmoothData, YsmoothData), axis=1)
                self.R0 = np.concatenate((self.R0, R0), axis=1)
                self.Y0 = np.concatenate((self.Y0, Y0), axis=1)
                self.Rgoodsc = np.concatenate((self.Rgoodsc, Rgoodsc), axis=1)
                self.Ygoodsc = np.concatenate((self.Ygoodsc, Ygoodsc), axis=1)
                self.R0sc = np.concatenate((self.R0sc, R0sc), axis=1)
                self.Y0sc = np.concatenate((self.Y0sc, Y0sc), axis=1)

            print(np.shape(self.Yscaled))
            print(np.shape(self.YsmoothData))
        
        self.getSNR(self.Yscaled, self.YsmoothData)
        self.getGoodComponents(snrTh, expFitTh, self.YsmoothData) # deprecated

        self.dYY = np.divide(self.Ygoodsc-self.Y0sc, self.Y0sc)
        self.dRR = np.divide(self.Rgoodsc-self.R0sc, self.R0sc)

        self.Ybl = np.transpose( np.transpose(self.Ygoodsc-self.Y0sc) + np.amin(self.Y0sc, axis=1) ) # bleach corrected but scaled like raw data
        self.Rbl = np.transpose( np.transpose(self.Rgoodsc-self.R0sc) + np.amin(self.R0sc, axis=1) )

        
        #self.dOO = np.divide(self.Ogoodsc-self.O0sc, self.O0sc)
        self.make_O_and_dOO()
        self.getGoodComponentsFull()
        # self.O = self.O[np.flatnonzero(self.oIsGood),:]
        # self.dOO = self.dOO[np.flatnonzero(self.oIsGood),:]
        # self.dYY = self.dYY[np.flatnonzero(self.oIsGood),:]
        # self.dRR = self.dRR[np.flatnonzero(self.oIsGood),:]
        
        # self.Opopt = self.popt
        self.calculate_I_and_dII()

        dataToCluster = self.dOO[np.flatnonzero(self.goodIds),:]
        self.computeCorr(dataToCluster)
        try:
            self.hierCluster(self.cc)
            self.plotClusters(dataToCluster) 
            self.computeICA(self.dOO, 6) #self.clustData
        except:
            self.cc = 0
            self.clustInd = 0
            self.clustData = 0
            self.S_ = 0
        
        self.saveSummary(outputFile)
        



# def expFun(x, a, b, c, d):
#     return a * np.exp(-b * x) + c*x + d

def polyFun(x, c, d):
    return c*x + d #a*x**3 + b*x**2 + 

def expFun(t, a, b, tau):
    return b + a*np.exp(- t/tau)

def dExpFun(t, a, tau):
    return  - a/tau*np.exp(- t/tau)

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



if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from scipy import sparse, optimize
    import numpy as np
    import pickle
    from scipy import io
    import smoothBehavior as beh

    baseFolder = '/Users/evan/Dropbox/_sandbox/sourceExtraction/good/_runAndFeed/190424_f3/' #feeding/190319_Trh_f1/' #_runAndFeed/190422_f1/' #
    behobj = beh.smoothData(baseFolder)
    behobj.getSmoothBeh()

    obj = scape(baseFolder)
    obj.postProcess('F.mat', 'post_fromYcc.mat')
    # obj.showFitSample(10)
    # obj.showFitSample(50)
    # obj.showFitSample(80)
    obj = scape(baseFolder)
    obj.postProcess('F_fromRed.mat', 'post_fromRcc.mat')
    # obj.showFitSample(10)
    # obj.showFitSample(50)
    # obj.showFitSample(80)
    # plt.plot(obj.Rpopt[:,2],obj.Ypopt[:,2],'.'); plt.show()
    



    # bnds=[[0,0,-100,-100], [10,0.2,0,100]] # for exp
    # bnds=[[-.000001,-.000001,-100,-100], [.000001,0.000001,0,100]] # for poly
    # bnds=[[-100,-100], [0,100]] # for poly
    # #self.makeExpF0(self.tot,bnds)
    # tFit = np.linspace(0, np.shape(self.tot)[0]-1, np.shape(self.tot)[0])
    # self.poptTot, pcov = optimize.curve_fit(expFun, tFit, self.tot, maxfev = 100000000, bounds=bnds )
    # self.F0tot = np.matrix(expFun(tFit, *self.poptTot))
    # print(np.shape(self.F0tot))
    # bnds[0][1]=0.7*self.poptTot[1]
    # bnds[1][1]=4.0*self.poptTot[1]
    # # self.normalizeRawF(self.Y)


    # # pOrd = 1
    # # self.makePolyF0(self.RsmoothData, pOrd)
    # # self.R0 = self.F0
    # # self.makePolyF0(self.YsmoothData, pOrd)
    # bnds=[[0,0,-100,-100], [10,0.2,0,100]]
    # self.makeExpF0(self.RsmoothData, bnds)
    # self.R0 = self.F0
    # self.makeExpF0(self.YsmoothData, bnds)
    # self.Y0 = self.F0

