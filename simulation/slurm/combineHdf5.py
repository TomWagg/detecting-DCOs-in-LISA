import numpy as np
import h5py as h5
import sys
import getopt

PATH_TO_FILES = "/n/holystore01/LABS/berger_lab/Lab/fbroekgaarden/DATA/all_dco_legacy_CEbug_fix/"

def combineMetallicities(zlist, dirname):
    #-- Input files to be combined

    print('combining the COMPASOutput files for the simulation in directory ', dirname)
    z0 = zlist[0]
    print(z0, 'z0')
    h5File_0 = h5.File(PATH_TO_FILES + str(dirname) + "/Z_" + str(z0) + "/STROOPWAFELcombined/COMPASOutput.h5", 'r')

    groups = list(h5File_0.keys())
    print(groups)

    #-- Output file
    h5File3 = h5.File(PATH_TO_FILES + str(dirname) + "/COMPASOutputCombined.h5", 'w')
    for group in groups:

        #-- Create a group in the output
        print('now at', group)
        if group=='pulsarEvolution':
            print(group)
        else:
            h5File3.create_group(group)

            datasets = list(h5File_0[group].keys())

            #-- Loop through datasets, combine and write
            for dataset in datasets:
                combined_dataset = h5File_0[group][dataset]
                for ind, z in enumerate(zlist[1:]):
                    print('at metallicity = ', z)

                    ###-- Concatenate datasets
                    h5File_z = h5.File(PATH_TO_FILES + str(dirname) +"/Z_" + str(z) + "/STROOPWAFELcombined/COMPASOutput.h5", 'r')
                    combined_dataset = np.concatenate([combined_dataset, h5File_z[group][dataset]])

                #-- Save datasets
                h5File3[group][dataset] = combined_dataset

# COMBINED HIGH Z AND OLD FOR TOM: 
zzlist = [0.0001, 0.00011, 0.00012, 0.00014, 0.00016, 0.00017, 0.00019, 0.00022, 0.00024, 0.00027, 0.0003, 0.00034, 0.00037, 0.00042, 0.00047, 0.00052, 0.00058, 0.00065, 0.00073, 0.00081, 0.0009, 0.00101, 0.00113, 0.00126, 0.0014, 0.00157, 0.00175, 0.00195, 0.00218, 0.00243, 0.00272, 0.00303, 0.00339, 0.00378, 0.00422, 0.00471, 0.00526, 0.00587, 0.00655, 0.00732, 0.00817, 0.00912, 0.01018, 0.01137, 0.01269, 0.01416, 0.01581, 0.01765, 0.01971, 0.022, 0.0244,  0.02705, 0.03, '0.0139', '0.014', '0.01411', '0.01421', '0.01432', '0.01443', '0.01453', '0.01464', '0.01475', '0.01486', '0.01497', '0.01508', '0.0152', '0.01531', '0.01542', '0.01554', '0.01565', '0.01577', '0.01589', '0.01601', '0.01613', '0.01625', '0.01637', '0.01649', '0.01661', '0.01674', '0.01686', '0.01699', '0.01711', '0.01724', '0.01737', '0.0175', '0.01763', '0.01776', '0.01789', '0.01803', '0.01816', '0.01829', '0.01843', '0.01857', '0.01871', '0.01885', '0.01899', '0.01913', '0.01927', '0.01941', '0.01956', '0.01971', '0.01985', '0.02']


def main(argv):
    dirname = "fiducial"
    try:
        opts, args = getopt.getopt(argv, "hd:", ["dirname="])
    except getopt.GetoptError:
        print('combineHdf5.py -d <dirname>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('combineHdf5.py -d <dirname>')
            sys.exit()
        elif opt in ("-d", "--dirname"):
            dirname = arg

    print(dirname)

    combineMetallicities(zlist=zzlist, dirname=dirname)
    print('completed')


if __name__ == "__main__":
    main(sys.argv[1:])
