import h5py as h5
import numpy as np
from os.path import isfile

import sys
sys.path.append("../src/")

from formation_channels import identify_formation_channels
from variations import variations

dt = np.dtype(float)
dtype = [("m_1", dt), ("m_2", dt), ("a_DCO", dt), ("e_DCO", dt),
         ("a_LISA", dt), ("e_LISA", dt), ("t_evol", dt), ("t_merge", dt),
         ("tau", dt), ("dist", dt), ("Z", dt), ("snr", dt), ("weight", dt),
         ("seed", dt), ("channel", np.dtype(int))]

runs = 50

floor_path = "/n/holystore01/LABS/berger_lab/Lab/fbroekgaarden/DATA/all_dco_legacy_CEbug_fix/"

# loop over all binary types and physics variations
for bt in ["BHBH", "BHNS", "NSNS"]:
    for v in range(len(variations)):
        full_data = None
        n_ten_year = np.array([], dtype=np.int)
        n_runs = 0
        missing = []

        # loop over each run
        for i in range(1, runs + 1):
            # check individual file exists
            fname = "{}_{}_{}.h5".format(bt, variations[v]["file"], i)
            if isfile(fname):
                n_runs += 1

                # open file and copy contents
                with h5.File(fname, "r") as f:
                    n_curr = f["simulation"].attrs["n_ten_year"].astype(np.int)
                    n_ten_year = np.concatenate((n_ten_year, n_curr))
                    if full_data is None:
                        full_data = f["simulation"][...].squeeze()
                    else:
                        add_data = f["simulation"][...].squeeze()
                        prev_len = len(full_data)
                        full_data.resize(prev_len + len(add_data))
                        full_data[prev_len:] = add_data
            else:
                missing.append(i)

        # as long as there is at least one file
        if len(missing) != runs:
            # let the user know which ones are currently missing
            print(v, len(missing), missing)

            # work out the formation channels
            model = variations[v]["file"]
            if model == "optimistic":
                model = "fiducial"
            floor_file = floor_path \
                + "{}/COMPASOutputReduced.h5".format(model)
            with h5.File(floor_file, "r") as floor:
                channels = identify_formation_channels(full_data["seed"],
                                                       floor)

            # write the rest of the files to a single file
            fname = "../data/{}_{}_all.h5".format(bt, variations[v]["file"])
            with h5.File(fname, "w") as file:
                file.create_dataset("simulation", (np.sum(n_ten_year),),
                                    dtype=dtype)
                for col, _ in dtype:
                    file["simulation"][col] = full_data[col] \
                        if col != "channel" else channels
                file["simulation"].attrs["n_ten_year"] = n_ten_year

        # otherwise make sure that the user knows no data is present
        else:
            print("No data found for {}".format(variations[v]["long"]))
