import numpy as np

def tophat_bandpass(nu_c, delta_nu, samples = 50):
    """Calculate a tophat bandpass for a given central frequency and width. 
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
    print(np.sum(weights * delta_nu))
    return (freqs, weights)

f,w = tophat_bandpass(100,30)




