#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 20:58:40 2019

@author: wudl
"""
import numpy as np

class Settings():
    """constants for calling
    """
    def __init__(self):
        """ Parameters
        
        - psi: radians
        - nu_dust_ref : GHz
        - Td : temperature of thermaldust in Kelvin
        - nu_c : central frequency of bandpass
        - delta_nu : width of bandpass
        - samples : number samples in bandpass
        - fwhm : (full width at half maximum) fwhm for smoothing
        note that, samples = delta_nu+1 for top_hat bandpass
        - D : optical aperture of telescope in unit meter
        - c : speed of light. unit : m/s
        """
        self.epsilon = 1
        self.psi = np.pi
        self.nu_dust_ref = 353  # GHz
        self.Td = 19.6  # K
        self.nu_sync_ref = 30 #GHz

        self.nu_c = 95   
        self.delta_nu = 10
        self.samples = 6 
        self.D = 0.72
        self.c = 299792458
        