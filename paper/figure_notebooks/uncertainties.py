import numpy as np
import h5py as h5
import astropy.units as u
import legwork

dco_types = ["BHBH", "BHNS", "NSNS"]
data_folder = "../../simulation/data/"
sim_folder = data_folder + "simulation_high_res_double/"

def get_ecc_uncertainty(t_obs=4*u.yr, source_threshold=7, harmonic_threshold=7, model="fiducial"):
    detectable_harmonics = {"BHBH": None, "BHNS": None, "NSNS": None}
    snr_uncertainty = {"BHBH": None, "BHNS": None, "NSNS": None}
    ecc_uncertainty = {"BHBH": None, "BHNS": None, "NSNS": None}
    max_harmonics = {"BHBH": None, "BHNS": None, "NSNS": None}

    for i, dco_type in enumerate(dco_types):
        with h5.File(sim_folder + "{}_{}_all.h5".format(dco_type, model), "r") as f:
            full_data = f["simulation"][...].squeeze()
            snr_mask = full_data["snr"] * np.sqrt(t_obs / (4 * u.yr)) > source_threshold
            data = full_data[snr_mask]

        sources = legwork.source.Source(m_1=data["m_1"] * u.Msun, m_2=data["m_2"] * u.Msun,
                                        dist=data["dist"] * u.kpc, a=data["a_LISA"] * u.AU,
                                        ecc=data["e_LISA"])

        detectable_harmonics[dco_type] = np.zeros(len(data)).astype(int)
        snr_uncertainty[dco_type] = np.zeros(len(data))
        max_harmonics[dco_type] = np.zeros(len(data)).astype(int)

        harmonics_required = sources.harmonics_required(sources.ecc)

        harmonic_groups = [(1, 10), (10, 100), (100, 1000), (1000, 10000)]
        for lower, upper in harmonic_groups:
            match = np.logical_and(harmonics_required > lower, harmonics_required <= upper)
            if match.any():
                snr_n_2 = legwork.snr.snr_ecc_stationary(m_c=sources.m_c[match],
                                                         f_orb=sources.f_orb[match],
                                                         ecc=sources.ecc[match],
                                                         dist=sources.dist[match],
                                                         t_obs=t_obs,
                                                         harmonics_required=upper,
                                                         interpolated_g=sources.g,
                                                         interpolated_sc=sources.sc,
                                                         ret_snr2_by_harmonic=True)

                # count harmonics above threshold
                detectable_harmonics[dco_type][match] = (snr_n_2**0.5 > harmonic_threshold).astype(int).sum(axis=1)

                max_harmonics[dco_type][match] = np.argmax(snr_n_2, axis=1) + 1

                # get the top two harmonics and sum them to get uncertainty
                top_snrs = np.sort(snr_n_2**(0.5), axis=1)[:, -2:]
                snr_uncertainty[dco_type][match] = 1 / top_snrs[:, 0] + 1 / top_snrs[:, -1]

        ecc_uncertainty[dco_type] = snr_uncertainty[dco_type]

    return detectable_harmonics, snr_uncertainty, ecc_uncertainty, max_harmonics

def sky_localisation(snr, fGW, L=2*u.AU):
    sigma_theta = 16.6 * (7 / snr) * (5e-4 * u.Hz / fGW) * (2 * u.AU / L) * u.deg
    return sigma_theta.to(u.deg)

def get_f_orb_uncertainty(snr, t_obs, f_orb):
    return (4 * np.sqrt(3) / np.pi / (snr * t_obs) / f_orb).decompose()

def get_f_orb_dot_uncertainty(snr, t_obs, f_orb_dot):
    return (6 * np.sqrt(5) / np.pi / (snr * t_obs**2) / f_orb_dot).decompose()

def get_Fprime_over_F(e):
    return e * (1256 + 1608 * e**2 + 111 * e**4) / (96 + 196 * e**2 - 255 * e**4 - 37 * e**6)

def get_m_c_uncertainty(f_orb, f_orb_dot, ecc, ecc_uncertainty, snr, t_obs):
    f_orb_uncertainty = get_f_orb_uncertainty(snr, t_obs, f_orb)
    f_orb_dot_uncertainty = get_f_orb_dot_uncertainty(snr, t_obs, f_orb_dot)

    return 11 / 5 * f_orb_uncertainty \
        + 3 / 5 * f_orb_dot_uncertainty \
        + 3 / 5 * get_Fprime_over_F(ecc) * ecc_uncertainty

def get_m_c_uncertainty_alt(f_orb, f_orb_dot, ecc, ecc_uncertainty, snr, t_obs):
    f_orb_uncertainty = get_f_orb_uncertainty(snr, t_obs, f_orb)
    f_orb_dot_uncertainty = get_f_orb_dot_uncertainty(snr, t_obs, f_orb_dot)

    return np.sqrt((11 / 5 * f_orb_uncertainty)**2 + (3 / 5 * f_orb_dot_uncertainty)**2 + (3 / 5 * get_Fprime_over_F(ecc) * ecc_uncertainty)**2)