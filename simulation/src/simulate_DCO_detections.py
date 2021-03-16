import h5py as h5
import numpy as np
import astropy.units as u
from scipy.integrate import quad
from scipy.special import lambertw
import getopt
import sys

from compas_processing import get_COMPAS_vars, mask_COMPAS_data
import legwork as lw

SNR_CUTOFF = 7
MW_SIZE = 100000


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
        opts, _ = getopt.getopt(sys.argv[1:], "hi:o:n:t:f", ["help",
                                                             "input=",
                                                             "output=",
                                                             "loops=",
                                                             "binary-type=",
                                                             "opt-flag"])
    except getopt.GetoptError as err:
        print(err)
        usage()
        sys.exit(2)

    # set default values
    input_filepath = 'COMPASOutput.h5'
    output_filepath = 'COMPASOutput_testing.h5'
    loops = 10
    binary_type = "BHNS"
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
        # mask only required DCOs
        dco_mask = mask_COMPAS_data(COMPAS_file, binary_type, (True, True,
                                                               pessimistic))

        # get all relevant variables
        compas_m_1, compas_m_2,\
            compas_Z, compas_a_DCO,\
            compas_e_DCO, compas_t_evol,\
            compas_weights,\
            compas_seeds = get_COMPAS_vars(COMPAS_file,
                                           "doubleCompactObjects",
                                           ["m_1, m_2, Metallicity1",
                                            "separationDCOFormation",
                                            "eccentricityDCOFormation",
                                            "tform", "weight", "seed"],
                                           dco_mask)

        # add units
        compas_m_1, compas_m_2 = compas_m_1 * u.Msun, compas_m_2 * u.Msun
        compas_a_DCO *= u.AU
        compas_t_evol *= u.Myr

    # work out metallicity bins
    compas_Z_unique = np.unique(compas_Z)
    inner_bins = np.array([compas_Z_unique[i]
                           + (compas_Z_unique[i+1] - compas_Z_unique[i])
                           / 2 for i in range(len(compas_Z_unique) - 1)])
    Z_bins = np.concatenate(([compas_Z_unique[0]], inner_bins,
                             [compas_Z_unique[-1]]))

    # create a random number generator
    rng = np.random.default_rng()

    # prep the temporary variable for parameters
    MAX_HIGH = 500
    dt = np.dtype(float)
    to_file = np.zeros(shape=(loops * MAX_HIGH,),
                       dtype=[("m_1", dt), ("m_2", dt), ("a_DCO", dt),
                              ("e_DCO", dt), ("a_LISA", dt),
                              ("e_LISA", dt), ("t_evol", dt),
                              ("t_merge", dt), ("tau", dt), ("dist", dt),
                              ("Z", dt), ("snr", dt), ("weight", dt),
                              ("seed", dt)])

    n_ten_year_list = np.zeros(loops)
    tot_ten = 0
    for milky_way in range(loops):
        # draw position parameters from Frankel Model
        tau, dist, Z_unbinned = simulate_mw(MW_SIZE)

        # work out COMPAS limits (and limit to Z=0.022)
        min_Z_compas = np.min(compas_Z_unique)
        max_Z_compas = np.max(compas_Z_unique[compas_Z_unique <= 0.022])

        # change metallicities above COMPAS limits to between solar and upper
        too_big = Z_unbinned > max_Z_compas
        Z_unbinned[too_big] = 10**(np.random.uniform(np.log10(0.01416),
                                                     np.log10(max_Z_compas),
                                                     len(Z_unbinned[too_big])))

        # change metallicities below COMPAS limits to lower limit
        too_small = Z_unbinned < max_Z_compas
        Z_unbinned[too_small] = min_Z_compas

        # sort by metallicity so everything matches up well
        Z_order = np.argsort(Z_unbinned)
        tau, dist, Z_unbinned = tau[Z_order], dist[Z_order],\
            Z_unbinned[Z_order]

        # bin the metallicities using Floor's bins
        h, _ = np.histogram(Z_unbinned, bins=Z_bins)

        # draw binaries for each metallicity bin, store indices
        binaries = np.zeros(MW_SIZE).astype(np.int)
        indices = np.arange(len(compas_m_1)).astype(np.int)
        total = 0
        for i in range(len(h)):
            if h[i] > 0:
                same_Z = compas_Z == compas_Z_unique[i]
                binaries[total:total + h[i]] = rng.choice(indices[same_Z],
                                                          h[i], replace=True)
                total += h[i]

        # TODO: remove this eventually
        if total != MW_SIZE:
            print(compas_Z_unique)
            print(Z_bins)
            print(np.sum(h), h)
            print(min_Z_compas, max_Z_compas)
            exit("PANIC: something funky is happening with the Z bins")

        # mask parameters for binaries
        m_1, m_2, a_DCO, e_DCO,\
            t_evol, w, Z, seed = compas_m_1[binaries], compas_m_2[binaries],\
            compas_a_DCO[binaries], compas_e_DCO[binaries],\
            compas_t_evol[binaries], compas_weights[binaries],\
            compas_Z[binaries], compas_seeds[binaries]

        # work out which binaries are still inspiralling
        t_merge = lw.evol.get_t_merge_ecc(ecc_i=e_DCO, a_i=a_DCO,
                                          m_1=m_1, m_2=m_2)
        insp = t_merge > (tau - t_evol)

        # trim out the merged binaries
        m_1, m_2, a_DCO, e_DCO,\
            t_evol, t_merge,\
            tau, dist, Z, w = m_1[insp], m_2[insp], a_DCO[insp], e_DCO[insp],\
            t_evol[insp], t_merge[insp], tau[insp], dist[insp], Z[insp],\
            w[insp], seed[insp]

        # evolve binaries to LISA
        e_LISA, a_LISA, f_orb_LISA = lw.evol.evol_ecc(ecc_i=e_DCO, a_i=a_DCO,
                                                      m_1=m_1, m_2=m_2,
                                                      t_evol=tau - t_evol,
                                                      n_step=2,
                                                      output_vars=["ecc", "a",
                                                                   "f_orb"])
        # we only care about the final state
        e_LISA = e_LISA[:, -1]
        a_LISA = a_LISA[:, -1]
        f_orb_LISA = f_orb_LISA[:, -1]

        sources = lw.source.Source(m_1=m_1, m_2=m_2, ecc=e_LISA, dist=dist,
                                   f_orb=f_orb_LISA)
        snr = sources.get_snr(verbose=True)

        ten_year = snr > (SNR_CUTOFF / np.sqrt(10 / 4))
        n_ten_year = len(snr[ten_year])
        n_ten_year_list[milky_way] = n_ten_year

        # store parameters in temporary variable
        to_file["m_1"][tot_ten:tot_ten + n_ten_year] = m_1[ten_year]
        to_file["m_2"][tot_ten:tot_ten + n_ten_year] = m_2[ten_year]
        to_file["a_DCO"][tot_ten:tot_ten + n_ten_year] = a_DCO[ten_year]
        to_file["e_DCO"][tot_ten:tot_ten + n_ten_year] = e_DCO[ten_year]
        to_file["a_LISA"][tot_ten:tot_ten + n_ten_year] = a_LISA[ten_year]
        to_file["e_LISA"][tot_ten:tot_ten + n_ten_year] = e_LISA[ten_year]
        to_file["t_evol"][tot_ten:tot_ten + n_ten_year] = t_evol[ten_year]
        to_file["tau"][tot_ten:tot_ten + n_ten_year] = tau[ten_year]
        to_file["dist"][tot_ten:tot_ten + n_ten_year] = dist[ten_year]
        to_file["Z"][tot_ten:tot_ten + n_ten_year] = Z[ten_year]
        to_file["snr"][tot_ten:tot_ten + n_ten_year] = snr[ten_year]
        to_file["weight"][tot_ten:tot_ten + n_ten_year] = w[ten_year]
        to_file["seed"][tot_ten:tot_ten + n_ten_year] = seed[ten_year]

        tot_ten += n_ten_year

    to_file = to_file[:tot_ten]

    # store all parameters in h5 file
    with h5.File(output_filepath, "w") as file:
        file.create_dataset("simulation", (tot_ten,),
                            dtype=[("m_1", dt), ("m_2", dt), ("a_DCO", dt),
                                   ("e_DCO", dt), ("a_LISA", dt),
                                   ("e_LISA", dt), ("t_evol", dt),
                                   ("t_merge", dt), ("tau", dt), ("dist", dt),
                                   ("Z", dt), ("snr", dt), ("weight", dt),
                                   ("seed", dt)])
        file["simulation"][...] = to_file
        file["simulation"].attrs["n_ten_year"] = n_ten_year


if __name__ == "__main__":
    main()
