import h5py as h5
import numpy as np
from os.path import isfile

import sys
sys.path.append("../src/")

from variations import variations

dt = np.dtype(float)
dtype = [("m_1", dt), ("m_2", dt), ("a_DCO", dt), ("e_DCO", dt),
         ("a_LISA", dt), ("e_LISA", dt), ("t_evol", dt), ("t_merge", dt),
         ("tau", dt), ("dist", dt), ("Z", dt), ("snr", dt), ("weight", dt),
         ("seed", dt)]

for bt in ["BHBH", "BHNS", "NSNS"]:
    for v in range(len(variations)):
        full_data = None
        n_ten_year = np.array([], dtype=np.int)
        n_runs = 0
        missing = []
        for i in range(1, 50 + 1):
            fname = "{}_{}_{}.h5".format(bt, variations[v]["file"], i)
            if isfile(fname):
                n_runs += 1
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

        if len(missing) != 50:
            print(v, len(missing), missing)
            fname = "../data/{}_{}_all.h5".format(bt, variations[v]["file"])
            with h5.File(fname, "w") as file:
                file.create_dataset("simulation", (np.sum(n_ten_year),),
                                    dtype=dtype)
                file["simulation"][...] = full_data
                file["simulation"].attrs["n_ten_year"] = n_ten_year
        else:
            print("No data found for {}".format(variations[v]["long"]))
