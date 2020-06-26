import numpy as np
import astropy.units as u
import astropy.constants as c
from scipy.special import jv

def chirp_mass(m1, m2):
    """
        Calculate the chirp mass of a binary
        
        Args:
            m1 --> [array_like, Msun] Primary mass
            m2 --> [array_like, Msun] Secondary mass
            
        Returns:
            Mc --> [array_like, Msun] Chirp masss
    """
    return (m1 * m2)**(3./5) / (m1 + m2)**(1./5)

def gw_frequency(a, m1, m2):
    """
        Calculate the gravitional wave frequency of a binary
        
        Args:
            a  --> [array_like, AU]   Semi-major axis
            m1 --> [array_like, Msun] Primary mass
            m2 --> [array_like, Msun] Secondary mass
            
        Returns:
            f  --> [array_like, Hz]   Gravitational wave frequency
    """
    f = 1 / np.pi * np.sqrt(c.G * (m1 + m2) / a**3)
    return f.to(u.Hz)

def F(e):
    """ 
        Calculate F(e) from Peters (1963) Eq.17
    
        Args:
            e    --> [array_like, unitless] eccentricity
            
        Returns:
            F(e) --> [array_like, unitless] function of eccentricity
    """
    return np.divide(1 + 73/24 * np.square(e) + 37/96 * np.power(e, 4), np.power(1 - np.square(e), 7/2))

def g(n, e):
    """ 
        Calculate g(n, e) from Peters (1963) Eq.20
        
        Args:
            n  --> [int]                  harmonic number
            e  --> [array_like, unitless] eccentricity
            
        Returns:
            g  --> [array_like, unitless] fraction of GW power in harmonic n        
    """
    ne = n * e
    return (n**4 / 32.0) \
        * (np.square(jv(n - 2, ne) - 2 * e * jv(n - 1, ne) + (2. / n) * jv(n, ne) + 2 * e * jv(n + 1, ne) - jv(n + 2, ne)) \
        + (1 - np.square(e)) * np.square(jv(n - 2, ne) - 2 * jv(n, ne) + jv(n + 2, ne)) \
        + (4 / (3.0 * n**2)) * np.square(jv(n, ne))
    )

def strain(Mc, D, f, n, e):
    """
        Calculate the strain for binaries
        
        Args:
            Mc  --> [array_like, Msun]     Chirp mass
            D   --> [array_like, kpc]      Luminosity distance to binary
            f   --> [array_like, Hz]       Gravitational wave frequency
            n   --> [int]                  Harmonic
            e   --> [array_like, unitless] Eccentricity
            
        Returns:
            hn  --> [array_like, unitless] strain in the nth harmonic
    """
    h = np.sqrt(2**(25/3) / 5) * (c.G * Mc)**(5./3) / (D * c.c**4) * ((f / 2) * np.pi)**(2./3) * np.sqrt(g(n, e)) / n
    return h.decompose()

def calc_snr_stationary(hn, f, Sn, Tobs=4*u.yr, nmax=100):
    """
        Calculate the signal-to-noise for (possibly) eccentric and stationary binaries
        
        Args:
            hn   --> [array_like, unitless] Strain at the first nmax harmonics
            f    --> [array_like, Hz]       Gravitational wave frequency
            Sn   --> [function]             Sensitivity curve
            Tobs --> [float, yr]            Observing time (default=LISA mission length)
            nmax --> [int]                  Maxmium number of harmonics to sum over
            
        Returns:
            snr  --> [array_like, unitless] Signal-to-noise ratio
    """
    noise_n = [np.sqrt(Sn(f / 2 * n)) for n in range(1, nmax)]
    snr_n = hn * np.sqrt(2 * Tobs) / (noise_n / u.Hz**(1./2))
    return np.sum(snr_n, axis=0).decompose()
                                      
def calc_snr(Mc, D, f, e, Sn, Tobs=4*u.yr, nmax=100):
    for n in range(1, 100):
        hcn = characteristic_strain(Mc, D, f, e, n)
    

def characteristic_strain(Mc, D, f, e, n):
    """ 
        Find the characteristic strain in the nth harmonic caused by a binary
        (See e.g. Breivik et al. 2019)
        
        Args:
            Mc  --> [array_like, Msun]     Chirp mass
            D   --> [array_like, kpc]      Distance
            f   --> [array_like, Hz]       Gravitational wave frequency (2 * orbital frequency)
            n   --> [int]                  Harmonic number
            
        Returns:
            hcn --> [array_like, unitless] Characteristic strain in nth harmonic
    """
    HC2_CONSTS = 2 * c.G**(5./3) / (3 * np.pi**(4./3) * c.c**3)
    hcn2 = HC2_CONSTS * np.power(Mc, 5./3) / np.power(D, 2) / np.power(f / 2. * n, 1./3) * \
           np.power(2. / n, 2./3) * (g(n, e) / F(e))
    return np.sqrt(hcn2).decompose()

def characteristic_strain_circ(Mc, D, f):
    return characteristic_strain(Mc, D, f, e=0, n=2)

def total_characteristic_strain(Mc, D, f, e, nmax=100):
    """ 
        Find the total characteristic strain caused by a binary in the 1st through nmax'th harmonics
        (See e.g. Breivik et al. 2019)
        
        Args:
            Mc   --> [array_like, Msun]     Chirp mass
            D    --> [array_like, kpc]      Distance
            f    --> [array_like, Hz]       Gravitational wave frequency (2 * orbital frequency)
            nmax --> [int]                  Maximum harmonic number
            
        Returns:
            hc   --> [array_like, unitless] Characteristic strain
    """
    return np.sum([characteristic_strain(Mc, D, f, e, n) for n in range(1, nmax)], axis=0)