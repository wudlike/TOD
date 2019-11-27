#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 15:13:53 2019

@author: wudl
"""
import healpy as hp
import numpy as np

class Data_t():
    
    def __init__(self,nside,ra,dev,Imap,Qmap=None,Umap=None):
        self.nside = nside
        self.ra = ra
        self.dev = dev
        
    def transfer_fun_K(self):
        K = 1
        return(K)
        
    def noise(self):
        n = 0
        return(n)
    
    def gain(self):
        g = 1
        return(g)
    
    def eff_antenna_area(self):
        A_e = 1
        return(A_e)
    
    def bandpass(self):
        F_v = 1
        return(F_v)
        
    def orientation_angle(self):
        phi = 0
        return(phi)
        
    def beam_shape(self):
        P = 1
        return(P)
    
    def cross_polar_leakage(self):
        epsilon = 0
        return(epsilon)
    
    def tod(self):
        index_map = hp.ang2pix(self.nside, theta=self.dev, phi=self.ra)
        d_t = hp.smoothing(Imap,fwhm=np.radians(1.))[index_map]
        return(d_t)

# import map    
data = hp.read_map(r'COM_CMB_IQU-commander_1024_R2.01_full.fits',h=True,field=(0,1,2))
Imap, Qmap, Umap = [data[0][i] for i in [0,1,2]]

# for single freq, e.g., 100GHz
#scan strategy
nside = 1024
ra = np.arange(0,2*np.pi,np.pi/5)
dev = np.ones(len(ra))*np.pi/2
data_t = Data_t.tod(nside,ra,dev,Imap)























