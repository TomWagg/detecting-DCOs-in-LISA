import numpy as np
import astropy.units as u
from scipy.integrate import quad
from scipy.special import lambertw

__all__ = ["simulate_mw"]


def draw_lookback_times(size, tm=12*u.Gyr, tsfr=6.8*u.Gyr):
    """Inverse CDF sampling of lookback times using Frankel+2018 Eq. 4

    Parameters
    ----------
    size : `int`
        How many samples to draw
    tsfr : `float`
        Model parameter fit in Frankel+2018
    tm : `float`, optional
        Age of Milky Way, by default 12*u.Gyr

    Returns
    -------
    tau: `float/array`
        Random lookback times
    """
    # draw a random uniform variable
    U = np.random.rand(size)

    # compute the normalising coefficient
    Nt = 1 / quad(lambda x: np.exp(-(tm.value - x) / tsfr.value),
                  0, tm.value)[0]

    # work out lookback times
    tau = tm + tsfr * np.log((Nt * tsfr.value + U - 1) / (Nt * tsfr.value))
    return tau


def draw_radii(size, t, alpha=0.3):
    """Inverse CDF sampling of lookback times using Frankel+2018 Eq. 5

    Parameters
    ----------
    size : `int`
        How many samples to draw
    alpha : `float`, optional
        Disc inside-out growth parameter, by default 0.3

    Returns
    -------
    R: `float/array`
        Random Galactocentric radius
    """
    # draw a random uniform variable
    U = np.random.rand(size)

    R = - 4 * u.kpc * (1 - alpha * (t / (8 * u.Gyr)))\
        * (lambertw((U - 1) / np.exp(1), k=-1).real + 1)
    return R


def draw_heights(size, zd=0.3*u.kpc):
    """Inverse CDF sampling of lookback times using McMillan 2011 Eq. 3

    Parameters
    ----------
    size : `int`
        How many samples to draw
    zd : `float`, optional
        Disc scale height, by default 0.3*u.kpc

    Returns
    -------
    z: `float/array`
        Random heights
    """
    # draw a random uniform variable
    U = np.random.rand(size)
    z = np.random.choice([-1, 1], len(U)) * zd * np.log(1 - U)
    return z

def get_metallicity(R, tau, tm=12*u.Gyr, Fm=-1, gradient=-0.075/u.kpc,
                    Rnow=8.7*u.kpc, gamma=0.3, zsun=0.0142):
    """Convert radius and time to metallicity using Frankel+2018 Eq.7 and
    Bertelli+1994 Eq.9

    Parameters
    ----------
    R : `float/array`
        Galactocentric radii
    tau : `float/array`
        Lookback times
    tm : `float`, optional
        Maximum lookback time, by default 12*u.Gyr
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

    Returns
    -------
    Z: `float/array`
        Metallicities corresponding to radii and times
    """
    FeH = Fm + gradient * R - (Fm + gradient * Rnow) * (1 - (tau / tm))**gamma
    return np.power(10, 0.977 * FeH + np.log10(zsun))


def distance_from_earth(R, z, theta, Rsun=8.2*u.kpc):
    """Convert radii, height and angle to a distance to the Earth using trig.

    Parameters
    ----------
    R : `float/array`
        Galactocentric radius
    z : `float/array`
        Height above Galactic plane
    theta : `float/array`
        Azimuthal angle in disc
    Rsun : `float`, optional
        Galactocentric radius of the sun, by default 8.2*u.kpc

    Returns
    -------
    D : `float/array`
        Distances
    """
    D = np.sqrt(z**2 + R**2 + Rsun**2 - 2 * Rsun * R * np.cos(theta))
    return D


def simulate_mw(n_binaries, tm=12 * u.Gyr, tsfr=6.8 * u.Gyr, alpha=0.3,
                zd=0.3 * u.kpc, Fm=-1, gradient=-0.075 / u.kpc,
                Rnow=8.7 * u.kpc, gamma=0.3, zsun=0.0142, Rsun=8.2 * u.kpc,
                ret_pos=False, lookback=True):
    """Draw a sample of birth times, distances and metallicities from a Milky
    Way model.

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

    tau = draw_lookback_times(n_binaries, tm=tm, tsfr=tsfr)
    R = draw_radii(n_binaries, tau, alpha=alpha)
    z = draw_heights(n_binaries, zd=zd)
    Z = get_metallicity(R, tau, tm=tm, Fm=Fm, gradient=gradient, Rnow=Rnow,
                        gamma=gamma, zsun=zsun)
    theta = 2 * np.pi * np.random.rand(n_binaries)
    D = distance_from_earth(R, z, theta)

    if lookback is False:
        tau = tm - tau

    if ret_pos:
        return tau, D, Z, (R, z, theta)
    else:
        return tau, D, Z
