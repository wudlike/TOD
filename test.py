#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 21:49:07 2019

@author: wudl
"""

#from components import CMB

#bb = CMB().BB(353,2.725)
#I_, Q_, U_ = CMB().cmb_map()
#import pysm
#a = pysm.common.convert_units('K_CMB','K_RJ',30.)

import numpy as np
import scipy.constants as constants

def B(nu, T):
    """Planck function. 
    :param nu: frequency in GHz at which to evaluate planck function.
    :type nu: float.
    :param T: temperature of black body. 
    :type T: float.
    :return: float -- black body brightness.
    """
    x = constants.h * nu * 1.e9 / constants.k / T
    return 2. * constants.h * (nu * 1.e9) ** 3 / constants.c ** 2 / np.expm1(x)

def dB(nu, T):
    """Differential planck function. 
    
    :param nu: frequency in GHz at which to evaluate differential planck function.
    :type nu: float.
    :param T: temperature of black body. 
    :type T: float.
    :return: float -- differential black body function. 
    """
    x = constants.h * nu * 1.e9 / constants.k / T
    return B(nu, T) / T * x * np.exp(x) / np.expm1(x)

def K_CMB2Jysr(nu): 
    """Kelvin_CMB to Janskies per steradian. Nu is in GHz.
    :param nu: frequency in GHz at which to calculate unit conversion.
    :type nu: float, numpy.ndarray
    :return: unit conversion coefficient - float.
    
    """
    return dB(nu, 2.7255) * 1.e26

def K_RJ2Jysr(nu):
    """Kelvin_RJ to Janskies per steradian. Nu is in GHz.
    :param nu: frequency in GHz at which to calculate unit conversion.
    :type nu: float, numpy.ndarray.                                                                                                                      
    :return: unit conversion coefficient - float. 
    
    """
    return  2. * (nu * 1.e9 / constants.c) ** 2 * constants.k * 1.e26




