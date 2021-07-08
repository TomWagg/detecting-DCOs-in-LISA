import sys
sys.path.append("../src/")
from variations import variations


def create_file(file):
    script_file = "scripts/combine_{}.sh".format(file)

    lines = []
    lines.append("#!/bin/bash\n")
    lines.append("#SBATCH -J combine_{}\n".format(file))
    lines.append("#SBATCH -n 6\n")
    lines.append("#SBATCH -N 1\n")
    lines.append("#SBATCH -p serial_requeue\n")
    lines.append("#SBATCH --mem 12000\n")
    lines.append("#SBATCH -t 00-3:00\n")
    lines.append("#SBATCH -o /n/home09/twagg/detecting-DCOs-in-LISA/simulation/slurm/logs/combine_{}_%A_%a.out\n".format(file))
    lines.append("#SBATCH -e /n/home09/twagg/detecting-DCOs-in-LISA/simulation/slurm/logs/combine_{}_%A_%a.err\n".format(file))
    lines.append("#SBATCH --mail-user thomas.wagg@cfa.harvard.edu\n")
    lines.append("#SBATCH --mail-type ALL\n")

    lines.append("module load Anaconda3/2020.11\n")
    lines.append("source activate lisa\n")

    main_line = 'python /n/home09/twagg/detecting-DCOs-in-LISA/simulation/slurm/combineHdf5.py --dirname {}'.format(file)
    lines.append(main_line)

    with open(script_file, "w") as f:
        f.writelines(lines)


for v in variations:
    if v["file"] in ["fiducial", "alpha0_1", "ccSNkick_100km_s", "ccSNkick_30km_s",
                     "maxNSmass2_0", "maxNSmass3_0", "noPISN", "rapid",
                     "wolf_rayet_multiplier_0_1", "wolf_rayet_multiplier_5"]:
        create_file(v["file"])