import numpy as np
import astropy.units as u
from scipy.integrate import quad
from scipy.special import lambertw
from scipy.stats import beta

__all__ = ["simulate_mw"]


def draw_lookback_times(size, component="thin_disc", tm=12*u.Gyr, tsfr=6.8*u.Gyr):
    """Inverse CDF sampling of lookback times. Thin and thick discs uses Frankel+2018 Eq.4,
    separated at 8 Gyr. The bulge matches the distribution in Fig.7 of Bovy+19 but accounts for sample's bias.

    Parameters
    ----------
    size : `int`
        How many samples to draw
    component: `str`
        Which component of the Milky Way
    tm : `float`, optional
        Age of Milky Way, by default 12*u.Gyr
    tsfr : `float`
        Model parameter fit in Frankel+2018

    Returns
    -------
    tau: `float/array`
        Random lookback times
    """
    if component == "thin_disc":
        U = np.random.rand(size)
        norm = 1 / quad(lambda x: np.exp(-(tm.value - x) / tsfr.value), 0, 8)[0]
        tau = tsfr * np.log((U * np.exp(tm / tsfr)) / (norm * tsfr.value) + 1)
    elif component == "thick_disc":
        U = np.random.rand(size)
        norm = 1 / quad(lambda x: np.exp(-(tm.value - x) / tsfr.value), 8, 12)[0]
        tau = tsfr * np.log((U * np.exp(tm / tsfr)) / (norm * tsfr.value) + np.exp(8 * u.Gyr / tsfr))
    elif component == "bulge":
        tau = beta.rvs(a=2, b=3, loc=6, scale=6, size=size) * u.Gyr
    return tau


def R_exp(t, alpha=0.3):
    return 4 * u.kpc * (1 - alpha * (t / (8 * u.Gyr)))


def draw_radii(size, R_0):
    """Inverse CDF sampling of galactocentric radii.

    Parameters
    ----------
    size : `int`
        How many samples to draw
    R_0 : `float`
        Scale length

    Returns
    -------
    R: `float/array`
        Random Galactocentric radius
    """
    U = np.random.rand(size)
    R = - R_0 * (lambertw((U - 1) / np.exp(1), k=-1).real + 1)
    return R


def draw_heights(size, z_d=0.3*u.kpc):
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
    U = np.random.rand(size)
    z = np.random.choice([-1, 1], size) * z_d * np.log(1 - U)
    return z


def get_metallicity(R, tau, tm=12*u.Gyr, Fm=-1, gradient=-0.075/u.kpc,
                    Rnow=8.7*u.kpc, gamma=0.3, zsun=0.0142):
    """Convert radius and time to metallicity using Frankel+2018 Eq.7 and Bertelli+1994 Eq.9

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


def simulate_mw(n_binaries, components=["thin_disc", "thick_disc", "bulge"],
                masses=[2.585e10, 2.585e10, 0.91e10]*u.Msun, tm=12 * u.Gyr, tsfr=6.8 * u.Gyr, alpha=0.3,
                Fm=-1, gradient=-0.075 / u.kpc, Rnow=8.7 * u.kpc, gamma=0.3, zsun=0.0142, Rsun=8.2 * u.kpc,
                ret_pos=False, lookback=True):
    """Draw a sample of birth times, distances and metallicities from a Milky
    Way model.

    Parameters
    ----------
    n_binaries : `int`
        Number of binaries to simulate
    components : `list of strings`
        Which components to include: any from 'thin_disc', 'thick_disc' and 'bulge'.
    masses : `list of floats`
        Corresponding masses to `components` and so must have the same length.
    tm : `float`, optional
        Maximum lookback time, by default 12*u.Gyr
    tsfr : `float`, optional
        Star formation timescale, by default 6.8*u.Gyr
    alpha : `float`, optional
        Disc inside-out growth parameter, by default 0.3
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
        Whether to return lookback time (uses birth time if False), by default True

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
    # work out the weight to assign to each component
    total_mass = np.sum(masses)
    mass_fractions = np.divide(masses, total_mass)

    # convert these weights to a number of binaries
    sizes = np.zeros(len(mass_fractions)).astype(int)
    for i in range(len(components) - 1):
        sizes[i] = np.round(mass_fractions[i] * n_binaries)
    sizes[-1] = n_binaries - np.sum(sizes)

    tau = [None for i in range(len(components))]
    R = [None for i in range(len(components))]
    z = [None for i in range(len(components))]

    # go through each component and get lookback time, radius and height
    for i, com in enumerate(components):
        tau[i] = draw_lookback_times(sizes[i], tm=tm, tsfr=tsfr, component=com)

        scale_length = R_exp(tau[i], alpha=alpha) if com == "thin_disc"\
            else 1/0.43 * u.kpc if com == "thick_disc"\
            else 1.5 * u.kpc

        R[i] = draw_radii(sizes[i], R_0=scale_length)

        scale_height = 0.3 * u.kpc if com == "thin_disc"\
            else 0.95 * u.kpc if com == "thick_disc"\
            else 1.5 * u.kpc

        z[i] = draw_heights(sizes[i], z_d=scale_height)

    # combine the samples
    tau = np.concatenate(tau)
    R = np.concatenate(R)
    z = np.concatenate(z)

    # compute the metallicity given the other values
    Z = get_metallicity(R, tau, tm=tm, Fm=Fm, gradient=gradient, Rnow=Rnow, gamma=gamma, zsun=zsun)

    # draw a random azimuthal angle uniformly
    theta = 2 * np.pi * np.random.rand(n_binaries)

    # convert parameters to a distance from the earth
    D = distance_from_earth(R, z, theta, Rsun=Rsun)

    # swap to time from start of MW if necessary
    if lookback is False:
        tau = tm - tau

    # return positions as well as rest if requested
    if ret_pos:
        return tau, D, Z, (R, z, theta)
    else:
        return tau, D, Z
