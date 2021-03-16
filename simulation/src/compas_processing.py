import numpy as np

__all__ = ["get_COMPAS_vars", "mask_COMPAS_data"]


def get_COMPAS_vars(compas_file, group, variables, mask=None):
    """Return a variable from the COMPAS output data

    Parameters
    ----------
    input_file : `hdf5 File`
        COMPAS file
    group : `str`
        Group within COMPAS file
    variables : `str` or `list`
        List of column names to access in group (or single column name)
    mask : `bool/array`
        Mask of binaries with same shape as each variable

    Returns
    -------
    var_list : `various`
        Single variable or list of variables (all masked)
    """
    var_list = None
    if isinstance(var_list, str):
        var_list = compas_file[group][variables][...].squeeze()
        if mask is not None:
            var_list = var_list[mask]
    else:
        if mask is not None:
            var_list = [compas_file[group][var][...].squeeze()[mask]
                        for var in variables]
        else:
            var_list = [compas_file[group][var][...].squeeze()
                        for var in variables]

    return var_list


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
