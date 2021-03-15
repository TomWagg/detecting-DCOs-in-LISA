import h5py as h5
import numpy as np
import astropy.units as u
from scipy.integrate import quad, odeint
from scipy.special import lambertw
import getopt
import sys

SNR_CUTOFF = 7
MW_SIZE = 100000


def get_COMPAS_variable(input_file, param):
    """Return a variable from the COMPAS output data

    Parameters
    ----------
    input_file : `hdf5 File`
        COMPAS file
    param : `tuple`
        Tuple of (name of variable column, hdf5 keyname)

    Returns
    -------
    variable : `various`
        The requested variable
    """
    xparam, fxparam = param
    return input_file[fxparam][xparam][...].squeeze()

def mask_COMPAS_data(input_file, DCO_type, flags):
    """Mask COMPAS data based on binary type and COMPAS flags

    Parameters
    ----------
    input_file : `hdf5 File`
        COMPAS file
    DCO_type : `{{ 'ALL', 'BHBH', 'BHNS', 'NSNS' }}`
        Double compact object type
    bool_mask : `tuple`
        Flags for masking (mask binaries not merging in a Hubble time,
                           mask binaries with RLOF secondary after CEE,
                           mask Pessimistic CE binaries)
    Returns
    -------
    mask : `bool/array`
        Mask that can be applied to COMPAS variables
    """
    hubble, RLOF, pessimistic = flags
    fDCO = input_file['doubleCompactObjects']

    # get the total number of binaries
    BINARIES = len(fDCO['stellarType1'][...].squeeze())

    # store the stellar type of both stars
    type1 = fDCO['stellarType1'][...].squeeze()
    type2 = fDCO['stellarType2'][...].squeeze()

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
        rlof_sec_post_CEE = fDCO['RLOFSecondaryAfterCEE'][...].squeeze()
        rlof_mask = np.logical_not(rlof_sec_post_CEE)
    else:
        rlof_mask = np.repeat(True, BINARIES)

    if pessimistic:
        opt_flag = fDCO['optimisticCEFlag'][...].squeeze()
        pessimistic_mask = np.logical_not(opt_flag)
    else:
        pessimistic_mask = np.repeat(True, BINARIES)

    # combine all masks
    mask = type_mask * hubble_mask * rlof_mask * pessimistic_mask

    return mask


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
