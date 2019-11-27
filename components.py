import pysm
import numpy as np
import healpy as hp

class CMB():
    '''Class defining attributes in TOD from cmb component
    
    The current possible attributes are:
    
    - 'cmb_I' : intensity of cmb from map -- numpy.ndarray or float.
    - 'cmb_Q' : Q map used -- numpy.ndarray or float
    - 'cmb_U' : U map used -- numpy.ndarray or float
    - unit: uk_CMB
    '''
    def __init__(self):
        pass
    
    def BB(self,nu,T):
        """Planck function
        """
        return pysm.common.B(nu,T)
    
    def cmb_map(self,config=None):
        """Return cmb_I, Q, U maps
        
        unit: uk_CMB
        """
        cmb_Map = hp.read_map(r'../data/psm/components/cmb/cmb_map.fits',
                       field=(0,1,2))
        cmb_I, cmb_Q, cmb_U = cmb_Map[:]
        return cmb_I, cmb_Q, cmb_U
    
    def dB(self,nu,T):
        """The derivative of planck function with respect to time
        
        i.e., dB(v,T)/dT
        """
        return pysm.common.dB(nu,T)

class ThermalDust():
    """The contribution from thermaldust to TOD
        
    unit: MJy/sr.
    """
    def __init__(self):
        pass
    
    def dust_ampl(self,config=None):
        """return ampl_I,Q,U maps from fits
        for partial sky: dust_ampl_I = dust_ampl[0][0:10]
        """
        
        dust_ampl = hp.read_map(r'../data/psm/components/thermaldust/thermaldust_ampl.fits',
                                field=(0,1,2))
        dust_ampl_I,dust_ampl_Q,dust_ampl_U = dust_ampl[:]
        return dust_ampl_I, dust_ampl_Q, dust_ampl_U
    
    def dust_ind(self):
        """Are I, Q, U the same?
        
        """
        dust_ind1 = hp.read_map(r'../data/psm/components/thermaldust/thermaldust_specind1.fits',
                                field=(0))
        return dust_ind1

class Synchrotron():
    """The contribution from synchrotron to TOD
    
    unit: MJy/sr
    """
    def __init__(self):
        pass
    
    def sync_ampl(self,config=None):
        sync_ampl = hp.read_map(r'../data/psm/components/synchrotron/synchrotron_ampl.fits',
                                field=(0,1,2))
        sync_ampl_I,sync_ampl_Q,sync_ampl_U = sync_ampl[:]
        return sync_ampl_I, sync_ampl_Q, sync_ampl_U
    
    def sync_ind(self):
        """unit: unitless"""
        sync_ind = hp.read_map(r'../data/psm/components/synchrotron/synchrotron_specind.fits',
                                field=(0))
        return sync_ind
        











