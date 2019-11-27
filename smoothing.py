import healpy as hp
import numpy as np
import settings
import pysm
import time
import argparse

parser = argparse.ArgumentParser(description='test running')
parser.add_argument('-td','--Add_td',default='True', choices=['True','False'],
                    help='smoothing thermaldust, default=True')
parser.add_argument('-sycn','--Add_sync',default='True', choices=['True','False'],
                    help='smoothing synchrotron, default=True')
parser.add_argument('-cmb','--Add_cmb',default='True', choices=['True','False'],
                    help='smoothing cmb, default=True')
args = parser.parse_args()
config = vars(args)

st = settings.Settings()

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
    
    def cmb_map(self,config=None):
        """Return cmb_I, Q, U maps
        
        unit: uk_CMB
        """
        cmb_Map = hp.read_map(r'../data/psm/components/cmb/cmb_map.fits',
                       field=(0,1,2))
        return cmb_Map

class Synchrotron():
    """The contribution from synchrotron to TOD
    
    unit: MJy/sr
    """
    def sync_ampl(self,config=None):
        sync_ampl = hp.read_map(r'../data/psm/components/synchrotron/synchrotron_ampl.fits',
                                field=(0,1,2))    
#        sync_ampl = hp.read_map(r'../data/psm/components/synchrotron/test_sync_ampl.fits',
#                                field=(0,1,2))          
        return sync_ampl
    
    def sync_ind(self):
        """unit: unitless"""
        sync_ind = hp.read_map(r'../data/psm/components/synchrotron/synchrotron_specind.fits',
                                field=(0))
#        sync_ind = hp.read_map(r'../data/psm/components/synchrotron/test_sync_ind.fits',
#                                field=(0))        
        return sync_ind
    
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
#        dust_ampl = hp.read_map(r'../data/psm/components/thermaldust/test_dust_ampl.fits',
#                                field=(0,1,2))        
        return dust_ampl
    
    def dust_ind(self):
        """Are I, Q, U the same?
        
        """
        dust_ind1 = hp.read_map(r'../data/psm/components/thermaldust/thermaldust_specind1.fits',
                                field=(0))
#        dust_ind1 = hp.read_map(r'../data/psm/components/thermaldust/test_dust_ind.fits',
#                                field=(0))
        return dust_ind1
    
def convert_units(unit1,unit2,nu):
    """
    Function to do unit conversions between Rayleigh-Jeans units, CMB
    units, and flux units.
    :param unit1: unit from which we are converting.
    :type unit1: str.
    :param unit2: unit to which we are converting.
    :type unit2: str.
    :param nu: frequency at which to calculate unit conversion.
    :type nu: float, np.ndarray.
    :returns: unit conversion coefficient - float, numpy.ndarray.
    """
    return pysm.common.convert_units(unit1,unit2,nu) 
    
def sync_map(nu,ampl,ind):
    """Return synchrotron maps"""
    return ampl*(nu/st.nu_sync_ref)**ind

def BB(nu,T):
    """Planck function
    - unit: W/(m^2.Hz)
    """
#    return pysm.common.B(nu,T)  # W/(m^2.Hz)
    return 10**20*pysm.common.B(nu,T)  # MJysr

def dB(nu,T):
    """The derivative of planck function with respect to time
    
    i.e., dB(v,T)/dT
    """
    return pysm.common.dB(nu,T)

def dust_map(nu,ampl,ind):
    """Reture dust maps"""
    return ampl*(nu/st.nu_dust_ref)**ind*BB(nu,st.Td)
    
def map_smoothing(map_in,fwhm):
    """Return smoothed map
    
    - fwhm : radians
    """
    return hp.smoothing(map_in,fwhm,verbose=False)

def FWHM(nu):
    """Return full width at half maximum for a given frequency
    in the unit of arcmin
    
    - nu : array_like with unit in GHz
    - D : optical aperture of telescope, 
    - fwhm : radians
    """
    lam = st.c/(nu*10**9) #unit: m
    fwhm = 1.22*lam/st.D #unit: rad 
#    fwhm_arcmin = fwhm*180/np.pi*60 # arcmin
    return fwhm

def bandpass2():
    """read bandpass from given normalized weights
    
    Freqs range from 75 GHz to 105 GHz
    """
    data_bp = np.loadtxt('../data/ali/bandpass/90_norm_sim.txt')
    freqs = data_bp[:,0]
    weights = data_bp[:,1]
    dnu = 1
    # normalize bandpass
    weights = data_bp[:,1]/sum(data_bp[:,1])
    return (freqs,weights),dnu

#cmb_I*convert_units('uK_CMB','MJysr',nu)

if __name__ == "__main__":
    t = time.time()
    # band provided by bandpass from ali_bandpass data
    print('Smoothing start')
    channels,dnu = bandpass2()
    band = channels[0]
    np.save('../data/psm/smoothed_comps/freq_range',band)
    """smoothing and saving maps
    
    Maps are saved as .npy format, one needs to use np.load() to load data
    """ 
    if config['Add_cmb'] == 'True':
        print("Smoothing cmb maps...")
        cmb_maps = np.array([map_smoothing(dB(nu,T=2.73)*CMB().cmb_map()+BB(nu,T=2.73),FWHM(nu)) for nu in band])
        print("Writing smoothed cmb maps to files...")
    if config['Add_sync'] == 'True':
        print("Smoothing synchrotron maps...")
        sy_ampl = Synchrotron().sync_ampl()  
        sy_ind = Synchrotron().sync_ind() 
        sync_maps = np.array([map_smoothing(sync_map(nu,sy_ampl,sy_ind),FWHM(nu)) for nu in band])
        print("Writing smoothed synchrotron maps to files...")
        np.save('../data/psm/smoothed_comps/synchrotron/sync_maps_smoothed',sync_maps)
    if config['Add_td'] == 'True':
        print("Smoothing thermaldust maps...")
        td_ampl = ThermalDust().dust_ampl()
        td_ind1 = ThermalDust().dust_ind()
        dust_maps = np.array([map_smoothing(dust_map(nu,td_ampl,td_ind1),FWHM(nu)) for nu in band])
        print("Writing smoothed thermaldust maps to files...")
        np.save('../data/psm/smoothed_comps/thermaldust/dust_maps_smoothed',dust_maps)
    e = time.time()
    print('Smoothing end')
    print('time cost :\n>>> %.4f seconds' %(e-t))
    
