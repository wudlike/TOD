import numpy as np
import healpy as hp

class ScanStratege:
    """
    scan stratege
    """
    def boresight():
        bs = np.loadtxt(r'../data/ali/BORESIGHT_TRACE/BS20201220C000.txt')
        return bs


#bs = ScanStratege.boresight()
#dec = bs[:,2]
#ra = bs[:,1]
#index_map = hp.ang2pix(nside=2048, theta=dec, phi=np.pi/2-ra)