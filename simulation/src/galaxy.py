import numpy as np
import astropy.units as u
from scipy.integrate import quad
from scipy.special import lambertw

__all__ = ["simulate_mw"]


def simulate_mw(n_binaries, tm=12 * u.Gyr, tsfr=6.8 * u.Gyr, alpha=0.3,
                zd=0.3 * u.kpc, Fm=-1, gradient=-0.075 / u.kpc,
                Rnow=8.7 * u.kpc, gamma=0.3, zsun=0.0142, Rsun=8.2 * u.kpc,
                ret_pos=False, lookback=True):
    """Summary

    Parameters
    ----------
    n_binaries : `int`
        Number of binaries to simulate
    tm : `float`, optional
        Maximum lookback time, by default 12*u.Gyr
    tsfr : `float`, optional
        Star formation timescale, by default 6.8*u.Gyr
    alpha : `float`, optional
        Disc inside-out growth parameter, by default 0.3
    zd : `float`, optional
        Disc scale height, by default 0.3*u.kpc
    Fm : `int`, optional
        Metallicity at centre of disc at tm, by default -1
    gradient : `float`, optional
        Metallicity gradient, by default -0.075/u.kpc
    Rnow : `float`, optional
        Radius at which present day metallicity is solar, by default 8.7*u.kpc
    gamma : `float`, optional
        Time dependence of chemical enrichment, by default 0.3
    zsun : `float`, optional
        Solar metallicity, by default 0.0142
    Rsun : `float`, optional
        Galactocentric radius of the sun, by default 8.2*u.kpc
    ret_pos : bool, optional
        Whether to return full positions or just distance, by default False
    lookback : bool, optional
        Whether to return lookback time (uses birth time if False), by default
        True

    Returns
    -------
    tau : `float/array`
        Lookback times
    D : `float/array`
        Distance to the Earth
    Z : `float/array`
        Metallicity
    pos : `tuple`
        Positions (R, z, theta), returned if ``ret_pos=True``
    """
    def random_lookback_times(U):
        # Inverse CDF sampling derived from Frankel+ 2018 Eq. 2
        Nt = 1 / quad(lambda x: np.exp(-(tm.value - x) / tsfr.value),
                      0, tm.value)[0]
        return tm + tsfr * np.log((Nt * tsfr.value + U - 1)
                                  / (Nt * tsfr.value))

    def random_birth_radii(U, t):
        # Inverse CDF sampling derived from Frankel+ 2018 Eq. 5
        return - 3 * u.kpc * (1 - alpha * (t / (8 * u.Gyr)))\
            * (lambertw((U-1)/np.exp(1), k=-1).real + 1)

    def random_heights(U):
        # Inverse CDF sampling derived from McMillan+ 2011 Eq. 3
        return np.random.choice([-1, 1], len(U)) * zd * np.log(1 - U)

    def metallicity(R, t):
        # Equation combining Bertelli+ 1994 Eq. 9 and Frankel+ 2018 Eq. 7
        return np.power(10, 0.977 * (Fm + gradient * R - (Fm + gradient * Rnow)
                                     * (1 - (t / tm))**(gamma))
                        + np.log10(zsun))

    def distance_from_earth(R, z, theta):
        return np.sqrt(z**2 + R**2 + Rsun**2 - 2 * Rsun * R * np.cos(theta))

    tau = random_lookback_times(np.random.rand(n_binaries))
    R = random_birth_radii(np.random.rand(n_binaries), tau)
    z = random_heights(np.random.rand(n_binaries))
    Z = metallicity(R, tau)
    theta = 2 * np.pi * np.random.rand(n_binaries)
    D = distance_from_earth(R, z, theta)

    if lookback is False:
        tau = tm - tau

    if ret_pos:
        return tau, D, Z, (R, z, theta)
    else:
        return tau, D, Z
