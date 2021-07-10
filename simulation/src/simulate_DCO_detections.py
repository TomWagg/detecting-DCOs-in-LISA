import h5py as h5
import numpy as np
import astropy.units as u
import getopt
import sys

from compas_processing import get_COMPAS_vars, mask_COMPAS_data
from galaxy import simulate_mw, simulate_simple_mw
import legwork as lw

SNR_CUTOFF = 7
MW_SIZE = 200000


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
        opts, _ = getopt.getopt(sys.argv[1:], "hi:o:n:t:fs", ["help",
                                                             "input=",
                                                             "output=",
                                                             "loops=",
                                                             "binary-type=",
                                                             "opt-flag",
                                                             "simple-mw"])
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
    use_simple_mw = False

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
        if option in ("-s", "--simple-mw"):
            use_simple_mw = True

    # open COMPAS file
    with h5.File(input_filepath, "r") as COMPAS_file:
        # mask only required DCOs
        dco_mask = mask_COMPAS_data(COMPAS_file, binary_type, (True, True, pessimistic))

        # get all relevant variables
        compas_m_1, compas_m_2,\
            compas_Z, compas_a_DCO,\
            compas_e_DCO, compas_t_evol,\
            compas_weights,\
            compas_seeds = get_COMPAS_vars(COMPAS_file,
                                           "doubleCompactObjects",
                                           ["M1", "M2", "Metallicity1",
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
    inner_bins = np.array([compas_Z_unique[i] + (compas_Z_unique[i+1] - compas_Z_unique[i]) / 2
                           for i in range(len(compas_Z_unique) - 1)])
    Z_bins = np.concatenate(([compas_Z_unique[0]], inner_bins,
                             [compas_Z_unique[-1]]))

    # create a random number generator
    rng = np.random.default_rng()

    # prep the temporary variable for parameters
    MAX_HIGH = 500
    dt = np.dtype(float)
    dtype = [("m_1", dt), ("m_2", dt), ("a_DCO", dt), ("e_DCO", dt),
             ("a_LISA", dt), ("e_LISA", dt), ("t_evol", dt), ("t_merge", dt),
             ("tau", dt), ("Z", dt), ("R", dt), ("z", dt), ("theta", dt), ("snr", dt), ("weight", dt),
             ("seed", dt)]
    to_file = np.zeros(shape=(loops * MAX_HIGH,), dtype=dtype)

    n_detect_list = np.zeros(loops)
    total_MW_weight = np.zeros(loops)
    tot_detect = 0
    for milky_way in range(loops):
        if not use_simple_mw:
            # draw parameters from Frankel Model
            tau, dist, Z_unbinned, pos = simulate_mw(MW_SIZE)
        else:
            # draw parameters from simple Milky Way (following Breivik+2020)
            tau, dist, Z_unbinned, pos = simulate_simple_mw(MW_SIZE)
        R, z, theta = pos

        # work out COMPAS limits (and limit to Z=0.022)
        min_Z_compas = np.min(compas_Z_unique)
        max_Z_compas = np.max(compas_Z_unique[compas_Z_unique <= 0.022])

        # change metallicities above COMPAS limits to between solar and upper
        too_big = Z_unbinned > max_Z_compas
        Z_unbinned[too_big] = 10**(np.random.uniform(np.log10(0.01416),
                                                     np.log10(max_Z_compas),
                                                     len(Z_unbinned[too_big])))

        # change metallicities below COMPAS limits to lower limit
        too_small = Z_unbinned < min_Z_compas
        Z_unbinned[too_small] = min_Z_compas

        # sort by metallicity so everything matches up well
        Z_order = np.argsort(Z_unbinned)
        tau, dist, Z_unbinned, R, z, theta = tau[Z_order], dist[Z_order], Z_unbinned[Z_order],\
            R[Z_order], z[Z_order], theta[Z_order]

        # bin the metallicities using Floor's bins
        h, _ = np.histogram(Z_unbinned, bins=Z_bins)

        # draw binaries for each metallicity bin, store indices
        binaries = np.zeros(MW_SIZE).astype(np.int)
        indices = np.arange(len(compas_m_1)).astype(np.int)
        total = 0
        for i in range(len(h)):
            if h[i] > 0:
                same_Z = compas_Z == compas_Z_unique[i]
                binaries[total:total + h[i]] = rng.choice(indices[same_Z], h[i], replace=True)
                total += h[i]

        # TODO: remove this eventually
        if total != MW_SIZE:
            print(compas_Z_unique)
            print(Z_bins)
            print(np.sum(h), h)
            print(min_Z_compas, max_Z_compas)
            exit("PANIC: something funky is happening with the Z bins")

        # mask parameters for binaries
        m_1, m_2, a_DCO, e_DCO, t_evol, w, Z, seed = compas_m_1[binaries], compas_m_2[binaries],\
            compas_a_DCO[binaries], compas_e_DCO[binaries],\
            compas_t_evol[binaries], compas_weights[binaries],\
            compas_Z[binaries], compas_seeds[binaries]

        # store the total weight of full population (for normalisation)
        total_MW_weight[milky_way] = np.sum(w)

        # work out which binaries are still inspiralling
        t_merge = lw.evol.get_t_merge_ecc(ecc_i=e_DCO, a_i=a_DCO, m_1=m_1, m_2=m_2)
        insp = t_merge > (tau - t_evol)

        # trim out the merged binaries
        m_1, m_2, a_DCO, e_DCO, t_evol, t_merge, tau, dist, Z, R, z, theta, w, seed = m_1[insp], m_2[insp],\
            a_DCO[insp], e_DCO[insp], t_evol[insp], t_merge[insp], tau[insp], dist[insp], Z[insp], R[insp],\
            z[insp], theta[insp], w[insp], seed[insp]

        # evolve binaries to LISA
        e_LISA, a_LISA, f_orb_LISA = lw.evol.evol_ecc(ecc_i=e_DCO, a_i=a_DCO, m_1=m_1, m_2=m_2,
                                                      t_evol=tau - t_evol, n_step=2,
                                                      output_vars=["ecc", "a", "f_orb"])
        # we only care about the final state
        e_LISA = e_LISA[:, -1]
        a_LISA = a_LISA[:, -1]
        f_orb_LISA = f_orb_LISA[:, -1]

        sources = lw.source.Source(m_1=m_1, m_2=m_2, ecc=e_LISA, dist=dist, f_orb=f_orb_LISA)
        snr = sources.get_snr(verbose=True)

        detectable = snr > SNR_CUTOFF
        n_detect = len(snr[detectable])
        n_detect_list[milky_way] = n_detect

        # store parameters in temporary variable
        to_file["m_1"][tot_detect:tot_detect + n_detect] = m_1[detectable]
        to_file["m_2"][tot_detect:tot_detect + n_detect] = m_2[detectable]
        to_file["a_DCO"][tot_detect:tot_detect + n_detect] = a_DCO[detectable]
        to_file["e_DCO"][tot_detect:tot_detect + n_detect] = e_DCO[detectable]
        to_file["a_LISA"][tot_detect:tot_detect + n_detect] = a_LISA[detectable]
        to_file["e_LISA"][tot_detect:tot_detect + n_detect] = e_LISA[detectable]
        to_file["t_evol"][tot_detect:tot_detect + n_detect] = t_evol[detectable]
        to_file["tau"][tot_detect:tot_detect + n_detect] = tau[detectable]
        to_file["R"][tot_detect:tot_detect + n_detect] = R[detectable]
        to_file["z"][tot_detect:tot_detect + n_detect] = z[detectable]
        to_file["theta"][tot_detect:tot_detect + n_detect] = theta[detectable]
        to_file["Z"][tot_detect:tot_detect + n_detect] = Z[detectable]
        to_file["snr"][tot_detect:tot_detect + n_detect] = snr[detectable]
        to_file["weight"][tot_detect:tot_detect + n_detect] = w[detectable]
        to_file["seed"][tot_detect:tot_detect + n_detect] = seed[detectable]

        tot_detect += n_detect

    to_file = to_file[:tot_detect]

    # store all parameters in h5 file
    with h5.File(output_filepath, "w") as file:
        file.create_dataset("simulation", (tot_detect,), dtype=dtype)
        file["simulation"][...] = to_file
        file["simulation"].attrs["n_detect"] = n_detect_list
        file["simulation"].attrs["total_MW_weight"] = total_MW_weight


if __name__ == "__main__":
    main()
