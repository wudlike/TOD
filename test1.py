import numpy as np 
import healpy as hp
import settings
import time

st = settings.Settings()
class ScanStratege:
    """
    scan stratege
    """
    def boresight(self):
        bs = np.loadtxt(r'../data/ali/BORESIGHT_TRACE/BS20201220C000.txt')
        dec = bs[:, 2]
        ra = bs[:, 1]
        index_map = hp.ang2pix(nside=2048, theta=dec, phi=np.pi/2-ra)
        return index_map, ang_phi

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

def tod1(index_map,phi):
    tot_maps_smoothed = cmb_maps[0]+sync_maps[0]+dust_maps[0] + (1-st.epsilon)/(1+st.epsilon)*((cmb_maps[1]+\
                    sync_maps[1]+dust_maps[1])[index_map]*np.cos(2*phi)+(cmb_maps[2]+\
                    sync_maps[2]+dust_maps[2])[index_map]*np.sin(2*phi))
    Tod = sum(np.multiply(np.mat(channels[1]).T*dnu,np.mat(map_combine_smoothed)))


cmb_maps = np.load(r'../data/psm/smoothed_comps/cmb/cmb_maps_smoothed.npy')
sync_maps = np.load(r'../data/psm/smoothed_comps/synchrotron/sync_maps_smoothed.npy')
dust_maps = np.load(r'../data/psm/smoothed_comps/thermaldust/dust_maps_smoothed.npy')


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

index_map, ang_phi = ScanStratege.boresight()

print('\nend reading\trunning...')
# bandpass from ali data
channels,dnu = bandpass2(data_bp_norm)
band = data_bp[:,0]
tod = tod(channels,dnu,band)[index_map]
hp.write_map('tod_sc.fits',tod,extra_header=([('name','TOD'),('nside','2048')]),overwrite=True)

t_end = time.time()
print('#'*8+' time cost '+'#'*8)
print(str((t_end-t_start)/60)+' : mins')
	