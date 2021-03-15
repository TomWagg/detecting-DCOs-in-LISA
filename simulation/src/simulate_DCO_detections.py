import csv
import h5py as h5
import numpy as np
import astropy.units as u
import astropy.constants as c
from scipy.integrate import quad, odeint
import scipy.interpolate as interpolate
from scipy.special import lambertw
import time
import getopt
import os
import sys

SNR_CUTOFF = 7
MW_SIZE = 100000

def interp2d_pairs(*args,**kwargs):
    """
        Same interface as interp2d but the returned interpolant will evaluate its inputs as pairs of values.
        Taken from https://stackoverflow.com/questions/35360756/scipy-interp2d-for-pairs-of-coordinates
    """
    
    # Internal function, that evaluates pairs of values, output has the same shape as input
    def interpolant(x, y, f):
        x, y = np.asarray(x), np.asarray(y)
        return (interpolate.dfitpack.bispeu(f.tck[0], f.tck[1], f.tck[2], f.tck[3], f.tck[4], x.ravel(), y.ravel())[0]).reshape(x.shape)
    
    # Wrapping the scipy interp2 function to call out interpolant instead
    return lambda x, y: interpolant(x, y, interpolate.interp2d(*args,**kwargs))

def interpolated_snr(Mc, D, e, forb, Tobs, interpolation):
    snr = np.sqrt(2**(25/3) / 5 * 2 * Tobs.to(u.s)) * (c.G * Mc)**(5/3) / c.c**4 * (np.pi * forb)**(2/3) / D * np.sqrt(interpolation(e, forb.value)) * u.Hz**(1/2)
    return snr.decompose()

def get_COMPAS_variable(input_file, param):
    """ 
        Return a variable from the COMPAS output data
        This function is adapted slightly from code provided by Floor Broekgaarden.
        
        Args:
            input_file   --> [h5]         Variable containing an h5 File
            param    --> [tuple]      Tuple of (name of variable column, hdf5 keyname)
            
        Returns:
            variable --> [array_like] The requested variable
            
        Example call:
            get_COMPAS_variable(input_file, ("M1", "doubleCompactObjects"))
    """
    xparam, fxparam = param
    return input_file[fxparam][xparam][...].squeeze()

def mask_COMPAS_data(input_file, DCO_type, flags):
    """ 
        Mask COMPAS data based on binary type and COMPAS flags
        This function is adapted slightly from code provided by Floor Broekgaarden.
        
        Args:
            input_file    --> [h5]                  Variable containing an h5 File
            DCO_type  --> [str]                 Type of binary: ['ALL', 'BHBH', 'BHNS', 'NSNS']
            bool_mask --> [tuple, boolean]      Flags: (mask binaries not merging in a Hubble time, 
                                                        mask binaries with RLOF secondary after CEE,
                                                        mask Pessimistic CE binaries)
                                               
        Returns:
            mask      --> [array_like, boolean] Mask that can be applied to COMPAS variables
    """
    hubble, RLOF, pessimistic = flags
    fDCO = input_file['doubleCompactObjects']
    
    # get the total number of binaries
    BINARIES = len(fDCO['stellarType1'][...].squeeze())
    
    # store the stellar type of both stars
    type1, type2 = fDCO['stellarType1'][...].squeeze(), fDCO['stellarType2'][...].squeeze()
    
    # create a mask on type (where BH=14 and NS=13)
    if DCO_type == "ALL":
        type_mask = np.repeat(True, BINARIES)
    elif DCO_type == "BHBH":
        type_mask = np.logical_and(type1 == 14, type2 == 14)
    elif DCO_type == "NSNS":
        type_mask = np.logical_and(type1 == 13, type2 == 13)
    elif DCO_type == "BHNS":
        type_mask = np.logical_or(np.logical_and(type1 == 14, type2 == 13),
                                  np.logical_and(type1 == 13, type2 == 14))
    else:
        print("Error: Invalid DCO_type")
        return
        
    # mask based on the Hubble time, RLOF and pessimistic flags
    if hubble:
        hubble_mask = fDCO['mergesInHubbleTimeFlag'][...].squeeze()
    else:
        hubble_mask = np.repeat(True, BINARIES)
        
    if RLOF:
        rlof_mask = np.logical_not(fDCO['RLOFSecondaryAfterCEE'][...].squeeze())
    else:
        rlof_mask = np.repeat(True, BINARIES)
        
    if pessimistic:
        pessimistic_mask = np.logical_not(fDCO['optimisticCEFlag'][...].squeeze())
    else:
        pessimistic_mask = np.repeat(True, BINARIES)
    
    # combine all masks
    mask = type_mask * hubble_mask * rlof_mask * pessimistic_mask
        
    return mask

def c0_peters(a0, e0):
    """
        Find constant c0 from Peters (1964) Eq. 5.11
        
        Args:
            a0 --> [array_like, AU]       initial semi-major axis
            e0 --> [array_like, unitless] initial eccentricity
            
        Returns:
            c0 --> [array_like, AU]       c0 constant from Eq. 5.11    
    """
    return (a0 * (1 - np.square(e0)) * np.power(e0,-12./19) * \
          np.power(1 + (121./304) * np.square(e0), -870./2299)).to(u.AU)

def beta_peters(m1, m2):
    """
        Find the constant beta from Peters (1964)
        
        Args:
            m1   --> [array_like, Msun]     Primary mass
            m2   --> [array_like, Msun]     Secondary mass
            
        Returns:
            beta --> [array_like, m^4 / s]  beta constant from Peters (1964)
    """
    beta = (64. / 5) * (c.G**3 * m1 * m2 * (m1 + m2)) / c.c**5
    return beta.to(u.m**4 / u.s)

def inspiral_time(e0, a0=None, m1=None, m2=None, c0=None, beta=None, small_e_tol=1e-2, large_e_tol=1 - 1e-2):
    """
        Calculate the inspiral time of a binary with Peters (1964) Eq. 5.14 and following approximations

        Args:
            e0          --> [list of floats] Initial eccentricity
            a0          --> [list of floats] Initial semi-major axis (required if c0 not provided)
            m1          --> [list of floats] Primary mass (required if beta not provided)
            m2          --> [list of floats] Secondary mass (required if beta not provided)
            c0          --> [list of floats] c0 constant (optional, leave as None to autocalculate)
            beta        --> [list of floats] Initial eccentricity
            small_e_tol --> [float]          Eccentricity below which to apply the small e approximation
            large_e_tol --> [float]          Eccentricity above which to apply the large e approximation

        Returns:
            t_inspiral  --> [list of floats] Time until merger
    """
    # check for valid input
    if a0 is None and c0 is None:
        raise ValueError("Must provide either a0 or c0")
    if (m1 is None or m2 is None) and beta is None:
        raise ValueError("Must provide either (m1, m2) or beta")

    # calculate c0 and beta if not provided
    if c0 is None:
        c0 = c0_peters(a0, e0)
    if beta is None:
        beta = beta_peters(m1, m2)

    def peters_5_14(e):
        """ Inspiral time from Peters Eq. 5.14 """
        return np.power(e, 29/19) * np.power(1 + (121/304)*e**2, 1181/2299) / np.power(1 - e**2, 3/2)

    if isinstance(e0, list) or isinstance(e0, np.ndarray):
        circular = e0 == 0.0
        small_e = np.logical_and(e0 < small_e_tol, e0 > 0.0)
        large_e = e0 > large_e_tol
        other_e = np.logical_not(np.logical_or(small_e, large_e))

        t_inspiral = np.zeros(len(e0)) * u.Gyr
        t_inspiral[circular] = a0[circular]**4 / (4 * beta[circular])
        t_inspiral[small_e] = c0[small_e]**4 / (4 * beta[small_e]) * e0[small_e]**(48/19)
        t_inspiral[large_e] = c0[large_e]**4 / (4 * beta[large_e]) * e0[large_e]**(48/19) \
                            * (768 / 425) * (1 - e0[large_e]**2)**(-1/2) * (1 + 121/304 * e0[large_e]**2)**(3480/2299)
        t_inspiral[other_e] = [((12 / 19) * c0[other_e][i]**4 / beta[other_e][i] * quad(peters_5_14, 0, e0[other_e][i])[0]).to(u.Gyr).value
                               for i in range(len(e0[other_e]))] * u.Gyr
    else:
        if e0 == 0.0:
            t_inspiral = a0**4 / (4 * beta)
        elif e0 < small_e_tol:
            t_inspiral = c0**4 / (4 * beta) * e0**(48/19)
        elif e0 > large_e_tol:
            t_inspiral = c0**4 / (4 * beta) * e0**(48/19) * (768 / 425) * (1 - e0**2)**(-1/2) * (1 + 121/304 * e0**2)**(3480/2299)
        else:
            t_inspiral = ((12 / 19) * c0**4 / beta * quad(peters_5_14, 0, e0)[0])
    return t_inspiral.to(u.Gyr)

def inspiral_time_quad(a0, e0, c0, beta):
    """
        Calculate the coalescence of a binary with Peters (1964) Eq. 5.14
        
        Args:
            a0      --> [array_like, AU]       Initial semi-major axis
            e0      --> [array_like, unitless] Initial eccentricity
            m1      --> [array_like, Msun]     Primary mass
            m2      --> [array_like, Msun]     Secondary mass
            
        Returns:
            t_inspr --> [array_like, Gyr]      Time from DCO formation to merger
    """
    
    def peters_5_14(e):
        """ Inspiral time from Peters Eq. 5.14 """
        return np.power(e, 29/19) * np.power(1 + (121/304)*e**2, 1181/2299) / np.power(1 - e**2, 3/2)
    
    if isinstance(e0, list) or isinstance(e0, np.ndarray):
        t_inspr = [((12 / 19) * c0[i]**4 / beta[i] * quad(peters_5_14, 0, e0[i])[0]).to(u.Gyr).value
               for i in range(len(e0))] * u.Gyr
    else:
        t_inspr = ((12 / 19) * c0**4 / beta * quad(peters_5_14, 0, e0)[0]).to(u.Gyr)
    return t_inspr

def dedt(e, t, beta, c0):
    """
        Rate of change of eccentricity from Peters (1964) Eq. 5.13
        
        Args:
            e    --> [array_like, unitless] eccentricity
            beta --> [array_like, m^4 / s]  beta constant from Peters (1964), see Eq. 5.9
            c0   --> [array_like, m]        c0 constant from Peters (1964), see Eq. 5.11
    
        Returns:
            dedt --> [array_like, 1 / s]    rate of change of eccentricity
    """
    return - (19.0 / 12.0) * (beta / np.power(c0, 4.0)) * (np.power(e, -29.0/19.0) \
            * np.power(1.0 - np.square(e), 3.0/2.0)) / np.power(1.0 + (121.0/304.0) * np.square(e), 1181.0/2299.0)

def e_to_a(e, c0):
    """
        Convert eccentricity to semi major axis using Peters (1964) Eq. 5.11
        
        Args:
            e  --> [array_like, unitless] eccentricity
            c0 --> [array_like, AU]       c0 constant from Peters (1964), see Eq. 5.11
            
        Returns:
            a  --> [array_like, AU]       semi-major axis
    """
    a = c0 / (1.0 - np.square(e)) * np.power(e, 12./19) * np.power(1 + (121./304 * np.square(e)), 870./2299)
    return a.to(u.AU)

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

def forb_from_a(a, m1, m2):
    """
        Calculate the orbital frequency of a list of binares using Kepler's third law
        (inverse of forb_from_a())

        Args:
            a    --> [list of floats] Semi-major axis
            m1   --> [list of floats] Primary mass
            m2   --> [list of floats] Secondary mass

        Returns:
            forb --> [list of floats] Orbital frequency
    """
    f = 1 / (2 * np.pi) * np.sqrt(c.G * (m1 + m2) / a**3)
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

def simulate_mw(Nbinaries, tm=12 * u.Gyr, tsfr=6.8 * u.Gyr, alpha=0.3, zd=0.3 * u.kpc, Fm=-1, 
                gradient=-0.075 / u.kpc, Rnow=8.7 * u.kpc, gamma=0.3, zsun=0.0142, Rsun=8.2 * u.kpc,
                ret_pos=False, lookback=True):
    """
        Simulate the history of the Milky Way using distributions from Frankel at al. 2018
        Most arguments refer to parameters from Frankel et al. 2018, see this paper for more details
    
        Args:
            Nbinaries --> Number of binaries to simulate
            tm        --> Maximum lookback time (Default: 12 Gyr)
            tsfr      --> Star formation timescale (Default: 6.8 Gyr)
            alpha     --> Disk inside-out growth parameter
            zd        --> Disk scale height (Default: 0.3 kpc)
            Fm        --> Metallicity of the gas at the centre of the disk at tm (Default: -1 dex)
            gradient  --> Metallicity gradient of the interstellar medium (Default: -0.075 dex/kpc)
            Rnow      --> Radius at which the present day birth metallicity is solar (Default: 8.7 kpc)
            gamma     --> Time dependence of chemical enrichment (Default: 0.3)
            zsun      --> Metallicity of the sun (Default: 0.0142)
            Rsun      --> Galactocentric radius of the sun (Default: 8.2 kpc)
            ret_pos   --> Boolean, whether to return positions as well as distances
            lookback  --> Boolean, whether to return lookback time (or birth time)

        Returns:
            tau       --> Lookback time (or birth time if lookback=False) in Gyr
            D         --> Distance to the Earth in kpc
            Z         --> Metallicity
            pos       --> Tuple of the radii and heights in kpc and the angles (only returned if ret_pos=True)
    """
    
    def random_lookback_times(U):
        # Inverse CDF sampling derived from Frankel+ 2018 Eq. 2
        Nt = 1 / quad(lambda x: np.exp(-(tm.value - x) / tsfr.value), 0, tm.value)[0]
        return tm + tsfr * np.log((Nt * tsfr.value + U - 1) / (Nt * tsfr.value))
    
    def random_birth_radii(U, t):
        # Inverse CDF sampling derived from Frankel+ 2018 Eq. 5
        return - 3 * u.kpc * (1 - alpha * (t / (8 * u.Gyr))) * (lambertw((U-1)/np.exp(1), k=-1).real + 1)
    
    def random_heights(U):
        # Inverse CDF sampling derived from McMillan+ 2011 Eq. 3
        return np.random.choice([-1, 1], len(U)) * zd * np.log(1 - U)
    
    def metallicity(R, t):
        # Equation combining Bertelli+ 1994 Eq. 9 and Frankel+ 2018 Eq. 7
        return np.power(10, 0.977 * (Fm + gradient * R - (Fm + gradient * Rnow) * (1 - (t / tm))**(gamma))
                        + np.log10(zsun))
    
    def distance_from_earth(R, z, theta):
        return np.sqrt(z**2 + R**2 + Rsun**2 - 2 * Rsun * R * np.cos(theta))

    tau = random_lookback_times(np.random.rand(Nbinaries))
    R = random_birth_radii(np.random.rand(Nbinaries), tau)
    z = random_heights(np.random.rand(Nbinaries))
    Z = metallicity(R, tau)
    theta = 2 * np.pi * np.random.rand(Nbinaries)
    D = distance_from_earth(R, z, theta)
    
    if lookback == False:
        tau = tm - tau
    
    if ret_pos:
        return tau, D, Z, (R, z, theta)
    else:
        return tau, D, Z

def usage():
    print("usage: python full_simulation.py [options]")
    print("\toptions:")
    print("\t\t-h, --help  : print usage instructions")
    print("\t\t-i, --input : path to COMPAS h5 input file")
    print("\t\t-n, --loops : number of simulations to run")
    print("\t\t-o, --output: path to output h5 file")

#####################################################################

def main():
    # get command line arguments and exit if error
    try:
        opts, _ = getopt.getopt(sys.argv[1:], "hi:o:n:t:f", ["help", "input=", "output=", "loops=", "binary-type=", "opt-flag"])
    except getopt.GetoptError as err:
        print(err)
        usage()
        sys.exit(2)

    # set default values
    input_filepath = '../data/COMPASOutput.h5'
    output_filepath = '../output/COMPASOutput_testing.h5'
    loops = 10
    binary_type = "BHNS"
    label = "-1"
    pessimistic = True

    # change defaults based on input
    for option, value in opts:
        if option in ("-h", "--help"):
            usage()
            return
        if option in ("-i", "--input"):
            input_filepath = value
        if option in ("-o", "--output"):
            output_filepath = value
        if option in ("-n" or "--loops"):
            loops = int(value)
        if option in ("-t", "--binary-type"):
            binary_type = value
        if option in ("-f", "--opt-flag"):
            pessimistic = False

    # open COMPAS file
    with h5.File(input_filepath, "r") as COMPAS_file:
        dco_mask = mask_COMPAS_data(COMPAS_file, binary_type, (True, True, pessimistic))

        compas_primary_mass = get_COMPAS_variable(COMPAS_file, ["M1", "doubleCompactObjects"])[dco_mask] * u.Msun
        compas_secondary_mass = get_COMPAS_variable(COMPAS_file, ["M2", "doubleCompactObjects"])[dco_mask] * u.Msun
        compas_chirp_mass = chirp_mass(compas_primary_mass, compas_secondary_mass)

        compas_Z = get_COMPAS_variable(COMPAS_file, ["Metallicity1", "doubleCompactObjects"])[dco_mask]
        compas_Z_unique = np.unique(compas_Z)

        inner_bins = np.array([compas_Z_unique[i] + (compas_Z_unique[i+1] - compas_Z_unique[i]) / 2 for i in range(len(compas_Z_unique) - 1)])
        Z_bins = np.concatenate(([compas_Z_unique[0]], inner_bins, [compas_Z_unique[-1]]))

        compas_aDCO = get_COMPAS_variable(COMPAS_file, ["separationDCOFormation", "doubleCompactObjects"])[dco_mask] * u.AU
        compas_eDCO = get_COMPAS_variable(COMPAS_file, ["eccentricityDCOFormation", "doubleCompactObjects"])[dco_mask]

        compas_t_evolve = get_COMPAS_variable(COMPAS_file, ["tform", "doubleCompactObjects"])[dco_mask] * u.Myr

        compas_weights = get_COMPAS_variable(COMPAS_file, ["weight", "doubleCompactObjects"])[dco_mask]

    f_range = np.logspace(-6, -2, 1000) * u.Hz
    e_range = np.linspace(0, 1, 1000)
    n_range = np.arange(1, 1000 + 1)

    snr_vals = np.load("../data/snr_interp.npy")
    snr_interpolation = interp2d_pairs(e_range, f_range, snr_vals)

    # create a random number generator
    rng = np.random.default_rng()

    # prep the temporary variable for parameters
    MAX_HIGH = 500
    dt = np.dtype(float)
    to_file = np.zeros(shape=(loops * MAX_HIGH,), dtype=[("m1", dt), ("m2", dt), ("aDCO", dt), ("eDCO", dt), \
                                                        ("aLISA", dt), ("eLISA", dt), ("te", dt), ("t_inspiral", dt), ("tau", dt), \
                                                        ("D", dt), ("Z", dt), ("snr", dt), ("weight", dt)])

    n_high_snr = np.zeros(loops)
    total_high_snr = 0
    for milky_way in range(loops):
        # draw position parameters from Frankel Model
        tau, D, Z_unbinned = simulate_mw(MW_SIZE)

        min_Z_compas = np.min(compas_Z_unique)
        max_Z_compas = np.max(compas_Z_unique[compas_Z_unique <= 0.022])
        Z_unbinned[Z_unbinned > max_Z_compas] = 10**(np.random.uniform(np.log10(0.01416), np.log10(max_Z_compas), len(Z_unbinned[Z_unbinned > max_Z_compas])))
        Z_unbinned[Z_unbinned < min_Z_compas] = min_Z_compas

        # sort by metallicity so everything matches up well
        Z_order = np.argsort(Z_unbinned)
        tau, D, Z_unbinned = tau[Z_order], D[Z_order], Z_unbinned[Z_order]

        # bin the metallicities using Floor's bins
        h, _ = np.histogram(Z_unbinned, bins=Z_bins)

        # draw correct number of binaries for each metallicity bin and store indices
        binaries = np.zeros(MW_SIZE).astype(np.int)
        indices = np.arange(len(compas_primary_mass)).astype(np.int)
        total = 0
        for i in range(len(h)):
            if h[i] > 0:
                binaries[total:total + h[i]] = rng.choice(indices[compas_Z == compas_Z_unique[i]], h[i], replace=True)
                total += h[i]

        if total != MW_SIZE:
            print(compas_Z_unique)
            print(Z_bins)
            print(np.sum(h), h)
            print(min_Z_compas, max_Z_compas)
            exit("PANIC: something funky is happening with the Z bins")

        # mask parameters for binaries
        m1, m2, aDCO, eDCO, t_evolve, w, Z = compas_primary_mass[binaries], compas_secondary_mass[binaries], compas_aDCO[binaries], compas_eDCO[binaries], \
                                             compas_t_evolve[binaries], compas_weights[binaries], compas_Z[binaries]

        c0 = c0_peters(aDCO, eDCO).to(u.m)
        beta = beta_peters(m1, m2).to(u.m**4 / u.s)
        
        # work out which binaries are still inspiralling
        t_inspiral = inspiral_time_quad(e0=eDCO, a0=aDCO, c0=c0, beta=beta)
        inspiraling = t_inspiral > (tau - t_evolve)

        m1, m2, beta, aDCO, eDCO, c0, t_evolve, t_inspiral, tau, D, Z, w = m1[inspiraling], m2[inspiraling], beta[inspiraling], aDCO[inspiraling], \
                                                                           eDCO[inspiraling], c0[inspiraling], t_evolve[inspiraling], t_inspiral[inspiraling], \
                                                                           tau[inspiraling], D[inspiraling], Z[inspiraling], w[inspiraling]

        # evolve eccentricity for inspiraling binaries
        eLISA = np.array([odeint(dedt, eDCO[i], [0, (tau[i] - t_evolve[i]).to(u.s).value], \
                        args=(beta[i].value, c0[i].value))[-1][0] for i in range(len(eDCO))])

        # convert to separation
        aLISA = e_to_a(eLISA, c0)

        # convert to frequency
        forb_LISA = forb_from_a(aLISA, m1, m2)
        forb_LISA[forb_LISA > 1 * u.Hz] = 0 * u.Hz

        # calculate the signal-to-noise ratio
        snr = interpolated_snr(chirp_mass(m1, m2), D, eLISA, forb_LISA, 4 * u.yr, interpolation=snr_interpolation)

        high_snr = snr > (SNR_CUTOFF / np.sqrt(10 / 4))
        high_snr_binaries = len(snr[high_snr])
        n_high_snr[milky_way] = high_snr_binaries

        # store parameters in temporary variable
        to_file["m1"][total_high_snr:total_high_snr + high_snr_binaries] = m1[high_snr]
        to_file["m2"][total_high_snr:total_high_snr + high_snr_binaries] = m2[high_snr]
        to_file["aDCO"][total_high_snr:total_high_snr + high_snr_binaries] = aDCO[high_snr]
        to_file["eDCO"][total_high_snr:total_high_snr + high_snr_binaries] = eDCO[high_snr]
        to_file["aLISA"][total_high_snr:total_high_snr + high_snr_binaries] = aLISA[high_snr]
        to_file["eLISA"][total_high_snr:total_high_snr + high_snr_binaries] = eLISA[high_snr]
        to_file["t_evolve"][total_high_snr:total_high_snr + high_snr_binaries] = t_evolve[high_snr]
        to_file["t_inspiral"][total_high_snr:total_high_snr + high_snr_binaries] = t_inspiral[high_snr]
        to_file["tau"][total_high_snr:total_high_snr + high_snr_binaries] = tau[high_snr]
        to_file["D"][total_high_snr:total_high_snr + high_snr_binaries] = D[high_snr]
        to_file["Z"][total_high_snr:total_high_snr + high_snr_binaries] = Z[high_snr]
        to_file["snr"][total_high_snr:total_high_snr + high_snr_binaries] = snr[high_snr]
        to_file["weight"][total_high_snr:total_high_snr + high_snr_binaries] = w[high_snr]

        total_high_snr += high_snr_binaries

    to_file = to_file[:total_high_snr]

    # store all parameters in h5 file
    with h5.File(output_filepath, "w") as file:
        file.create_dataset("simulation", (total_high_snr,), dtype=[("m1", dt), ("m2", dt), ("aDCO", dt), ("eDCO", dt), \
                                                                   ("aLISA", dt), ("eLISA", dt), ("t_evolve", dt), ("t_inspiral", dt), \
                                                                   ("tau", dt), ("D", dt), ("Z", dt), ("snr", dt), ("weight", dt)])
        file["simulation"][...] = to_file
        file["simulation"].attrs["n_high_snr"] = n_high_snr

if __name__ == "__main__":
    main()
