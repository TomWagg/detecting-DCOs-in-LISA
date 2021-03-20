import numpy as np
from collections import defaultdict

import sys
sys.path.append("../src/")
from compas_processing import get_COMPAS_vars


def find_mass_transfer_moments(seeds):
    """Create an array with each moment of mass transfer labelled. This is done
    such that the first moment is labelled '1' at its index and second as '2'
    etc.

    Parameters
    ----------
    seeds : `int/array`
        List of seeds that correspond to a binary and each instance
        corresponds to a moment of mass transfer.

    Returns
    -------
    moments : `int/array`
        List same shape as ``seeds`` with each index labelled with which moment
        of mass transfer that the binary is going through.
    """
    # start a counter that defaults to 0
    counter = defaultdict(int)
    moments = np.zeros_like(seeds).astype(int)

    # loop over every seed
    for i in range(len(seeds)):
        # increment seed's count
        counter[seeds[i]] += 1

        # label index with mass transfer moment (1 is first)
        moments[i] = counter[seeds[i]]
    return moments


def identify_formation_channels(seeds, file):
    """Identify the formation channel that produced each seed. We consider 5
    main channels: classic, only stable, single core CEE, double core CEE and
    other. We define the channels as follows (the numbers are what is put in
    the ``channels`` output):

    Classic (1) -- In the first mass transfer, the primary overflows its
    Roche lobe stably and has a clear core-envelope structure (HG - TPAGB)
    whilst the secondary star is still on the main sequence. In the second mass
    transfer, the secondary overflows its Roche Lobe before stripping occurs
    and leads to a common envelope event.

    Only stable (2) -- As Classic, but the second mass transfer is stable.

    Single core CEE (3) -- In the first mass transfer, the primary overflows
    its Roche lobe unstably (leading to a CEE) and has a clear core-envelope
    structure (FGB - TPAGB) whilst the secondary star is still on the main
    sequence.

    Double core CEE (4) -- In the first mass transfer, the primary overflows
    its Roche lobe unstably (leading to a CEE) and both the primary and
    secondary have a clear core-envelope structure (FGB - TPAGB).

    Other (5) -- Anything else

    Parameters
    ----------
    seeds : `int/array`
        List of seeds that each correspond to a binary (this should be a subset
        of the seeds in the COMPAS output file)
    file : `hdf5 File`
        An open hdf5 file (returned by h5py.File) with the COMPAS output

    Returns
    -------
    channels : `int/array`
        List of channels through with each binary formed
    """
    all_rlof_seeds = get_COMPAS_vars(file, "RLOF", "randomSeed")
    rlof_mask = np.isin(all_rlof_seeds, seeds)
    rlof_seeds = all_rlof_seeds[rlof_mask]

    moments = find_mass_transfer_moments(rlof_seeds)

    rlof_primary, rlof_secondary, cee_flag,\
        stellar_type_1, stellar_type_2 = get_COMPAS_vars(file,
                                                         "RLOF",
                                                         ["flagRLOF1",
                                                          "flagRLOF2",
                                                          "flagCEE",
                                                          "type1Prev",
                                                          "type2Prev"],
                                                         rlof_mask)

    # CLASSIC channel
    # 1st transfer, stable RLOF from primary (post-MS, unstripped) onto MS
    classic_or_OS_MT1 = np.logical_and.reduce((moments == 1,
                                               rlof_primary,
                                               np.logical_not(cee_flag),
                                               stellar_type_1 > 1,
                                               stellar_type_1 < 7,
                                               stellar_type_2 <= 1))
    # 2nd transfer, unstripped secondary RLOF into CE
    classic_MT2 = np.logical_and.reduce((moments == 2,
                                         rlof_secondary,
                                         cee_flag,
                                         stellar_type_2 < 7))
    classic_seeds = np.intersect1d(rlof_seeds[classic_or_OS_MT1],
                                   rlof_seeds[classic_MT2])
    classic_mask = np.isin(seeds, classic_seeds)

    # ONLY STABLE channel
    # 1st transfer as classic, 2nd transfer unstripped secondary stable RLOF
    only_stable_MT2 = np.logical_and.reduce((moments == 2,
                                             rlof_secondary,
                                             np.logical_not(cee_flag),
                                             stellar_type_2 < 7))
    only_stable_seeds = np.intersect1d(rlof_seeds[classic_or_OS_MT1],
                                       rlof_seeds[only_stable_MT2])
    only_stable_mask = np.isin(seeds, only_stable_seeds)

    # SINGLE CORE CEE channel
    # 1st transfer unstable, primary giant branch onto MS secondary
    single_core = np.logical_and.reduce((moments == 1,
                                         rlof_primary,
                                         cee_flag,
                                         stellar_type_1 > 2,
                                         stellar_type_1 < 7,
                                         stellar_type_2 <= 1))
    single_core_mask = np.isin(seeds, rlof_seeds[single_core])

    # DOUBLE CORE CEE channel
    # 1st transfer unstable, primary giant branch onto giant branch secondary
    double_core = np.logical_and.reduce((moments == 1,
                                         rlof_primary,
                                         cee_flag,
                                         stellar_type_1 > 2,
                                         stellar_type_1 < 7,
                                         stellar_type_2 > 2,
                                         stellar_type_2 < 7))
    double_core_mask = np.isin(seeds, rlof_seeds[double_core])

    channels = np.zeros_like(seeds).astype(int)
    channels[classic_mask] = 1
    channels[only_stable_mask] = 2
    channels[single_core_mask] = 3
    channels[double_core_mask] = 4

    return channels
