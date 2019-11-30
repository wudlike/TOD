import healpy as hp
import numpy as np
import settings
import pysm
import time

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

def BB(nu,T):
    """Planck function
       return black body brightness
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

if __name__ == "__main__":
    t = time.time()
    # band provided by bandpass from ali_bandpass data
    print('Smoothing start')
    channels,dnu = bandpass2()
    # band = channels[0]
    band = [90,91]
    np.save('../data/psm/smoothed_comps/freq_range',band)
    """smoothing and saving maps
    
    Maps are saved as .npy format, one needs to use np.load() to load data
    """ 
    print("Smoothing cmb maps...")
    cmb_map_in = CMB().cmb_map()
    cmb_maps = np.array([map_smoothing(dB(nu,T=2.73)*cmb_map_in*convert_units('uK_CMB','MJysr',nu)+BB(nu,T=2.73),FWHM(nu)) for nu in band])
    print("Writing smoothed cmb maps to files...")
    np.save('../data/psm/smoothed_comps/cmb/cmb_maps_smoothed_test',cmb_maps)