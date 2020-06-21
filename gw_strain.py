import numpy as np
import astropy.units as u

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
    hcn2 = HC2_CONSTS * np.power(Mc, 5./3) / np.power(D, 2) / np.power(f / 2 * n, 1./3) * \
           np.power(2 / n, 2./3) * (g(n, e) / F(e))
    return np.sqrt(hcn2).decompose()

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