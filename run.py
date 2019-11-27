from components import CMB,ThermalDust,Synchrotron
from scan_stratege import ScanStratege
import healpy as hp
import numpy as np
import pysm
import time
import settings

st = settings.Settings()

def display_time(func):
    """show running time
    """
    def wrapper(*args):
        t_start = time.time()
        TOD = func(*args)
        t_end = time.time()
        print('#'*6+' time cost '+'#'*6)
        print(str((t_end-t_start)/60)+' : mins')
        return TOD
    return wrapper

def dust_map(nu,ind,ampl):
    """Reture dust maps"""
    return ampl*(nu/st.nu_dust_ref)**ind*(CMB().BB(nu,st.Td)/CMB().BB(st.nu_dust_ref,st.Td))**2

def sync_map(nu,ind,ampl):
    """Return synchrotron maps"""
    return ampl*(nu/st.nu_sync_ref)**ind

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
    
def tophat_bandpass(nu_c, delta_nu, samples = 50):
    """Calculate a normalized tophat bandpass for a given central frequency and width. 
    This will return a tuple containing (frequencies, weights). 
    :param nu_c: central frequency of bandpass.
    :type nu_c: float.
    :param delta_nu: width of bandpass.
    :type delta_nu: float.
    :param samples: number samples in bandpass; the more samples the more accurate the result.
    :type samples: int.
    :return: tuple - (frequencies, weights)
    """
    freqs = np.linspace(nu_c - delta_nu / 2., nu_c + delta_nu / 2., samples)
    weights = np.ones_like(freqs) / (freqs.size * delta_nu / samples)
    return (freqs, weights)

def bandpass1(nu_c,delta_nu,samples):
    """I_Mjysr(nu)*weights*dnu
    """
    # tophat_bandpass
    dnu = delta_nu/(samples-1)
    band_pass = tophat_bandpass(nu_c = nu_c, delta_nu = delta_nu, samples = samples)
    return band_pass,dnu

def bandpass2(data):
    """read bandpass from given data"""
    freqs = data[:,0]
    weights = data[:,1]
    dnu = 1
    return (freqs,weights),dnu

def map_smoothing(nu,map_in,fwhm):
    """Return smoothed map
    
    - fwhm : radians
    """
    return hp.smoothing(map_in,fwhm,verbose=False)

#@display_time
def tod(channels,dnu,band):
    """Return time ordered data(TOD)
    
    -map_combine_smoothed : smoothing the whole combined component maps(cmb+dust+sync+...)
                   smooting([I+(1-epsilon)/(1+epsilon)(Q...U)])
    -TOD : ndarray
    -st.psi : radians
    """
    map_combine_smoothed = [map_smoothing(nu,cmb_I*convert_units('uK_CMB','MJysr',nu)\
                                          +dust_map(nu,dust_ind1,dust_ampl_I)+sync_map(nu,sync_ind,sync_ampl_I)+\
                            (1-st.epsilon)/(1+st.epsilon)*((cmb_Q*convert_units('uK_CMB','MJysr',nu)+\
                                            dust_map(nu,dust_ind1,dust_ampl_Q)+\
                                            sync_map(nu,sync_ind,sync_ampl_Q))*np.cos(2*st.psi)+\
                            (cmb_U*convert_units('uK_CMB','MJysr',nu)+dust_map(nu,dust_ind1,dust_ampl_U)+\
                             sync_map(nu,sync_ind,sync_ampl_U))*\
                            np.sin(2*st.psi)),FWHM(nu)) for nu in band]
    Tod = sum(np.multiply(np.mat(channels[1]).T*dnu,np.mat(map_combine_smoothed)))
#    print(Tod)
    return np.array(Tod)[0]

def conv(channels,dnu,band):
    """save I,Q,U maps after smoothing and convolved with bandpass
    """
    map_I_smoothed = [map_smoothing(nu,dust_map(nu,dust_ind1,dust_ampl_I),FWHM(nu)) for nu in band]
    I_conv = sum(np.multiply(np.mat(channels[1]).T*dnu,np.mat(map_I_smoothed)))
    map_Q_smoothed = [map_smoothing(nu,dust_map(nu,dust_ind1,dust_ampl_I),FWHM(nu)) for nu in band]
    Q_conv = sum(np.multiply(np.mat(channels[1]).T*dnu,np.mat(map_Q_smoothed)))
    map_U_smoothed = [map_smoothing(nu,dust_map(nu,dust_ind1,dust_ampl_I),FWHM(nu)) for nu in band]
    U_conv = sum(np.multiply(np.mat(channels[1]).T*dnu,np.mat(map_U_smoothed)))
    return (np.array(I_conv),np.array(Q_conv),np.array(U_conv))

def save_file(file):
    """save tod as fits format
    """
    return hp.write_map('I_conv.fits',file,extra_header=([('name','TOD'),('nside','2048')]),overwrite=True)
#    return hp.write_map(str(file)+'.fits',file,extra_header=([('name',str(file)),('nside','2048')]),overwrite=True)


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

t_start = time.time()
print('\nreading data ...\n')
dust_ampl_I, dust_ampl_Q, dust_ampl_U = ThermalDust().dust_ampl()
dust_ind1 = ThermalDust().dust_ind()
sync_ampl_I, sync_ampl_Q, sync_ampl_U = Synchrotron().sync_ampl()
sync_ind = Synchrotron().sync_ind()
cmb_I = CMB().cmb_map()[0]
cmb_Q = CMB().cmb_map()[1]
cmb_U = CMB().cmb_map()[2]
data_bp = np.loadtxt('../data/ali/bandpass/90_norm_sim.txt')
# normalize bandpass
data_bp = data_bp[:,:]
data_bp[:,1] = data_bp[:,1]/sum(data_bp[:,1])
data_bp_norm = data_bp

bs = ScanStratege.boresight()
dec = bs[:,2]
ra = bs[:,1]
index_map = hp.ang2pix(nside=1, theta=dec, phi=np.pi/2-ra)

#nside = 1
#dust_ampl_I = np.arange(12*nside**2)  # nside
#dust_ampl_Q = dust_ampl_I
#dust_ampl_U = dust_ampl_I
#dust_ind1 = np.arange(12*nside**2)*0.1
#
#sync_ampl_I = np.arange(12*nside**2)
#sync_ampl_Q = sync_ampl_I
#sync_ampl_U = sync_ampl_I
#sync_ind = np.arange(12*nside**2)*0.1
#cmb_I = np.arange(12*nside**2)
#cmb_Q = cmb_I
#cmb_U = cmb_I
print('\nend reading\trunning...')
# tophat_bandpass
#channels,dnu = bandpass1(st.nu_c,st.delta_nu,st.samples)  # return sum(delta_nu*weight)
#band = np.linspace(st.nu_c-st.delta_nu/2,st.nu_c+st.delta_nu/2,st.samples)
# bandpass from ali data
channels,dnu = bandpass2(data_bp_norm)
band = data_bp[:,0]
tod = tod(channels,dnu,band)[index_map]
hp.write_map('tod_sc.fits',tod,extra_header=([('name','TOD'),('nside','2048')]),overwrite=True)
#map_conv = conv(channels,dnu,band)
#hp.write_map('I_conv.fits',map_conv[0],extra_header=([('name','I_map_conv'),('nside','2048')]),overwrite=True)
#hp.write_map('Q_conv.fits',map_conv[1],extra_header=([('name','Q_map_conv'),('nside','2048')]),overwrite=True)
#hp.write_map('U_conv.fits',map_conv[2],extra_header=([('name','U_map_conv'),('nside','2048')]),overwrite=True)
t_end = time.time()
print('#'*8+' time cost '+'#'*8)
print(str((t_end-t_start)/60)+' : mins')












