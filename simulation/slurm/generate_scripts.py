import sys
sys.path.append("../src/")
from variations import variations


def create_file(bt, file, simple_mw=False):
    script_file = "scripts/{}_{}.sh".format(bt, file)

    lines = []
    lines.append("#!/bin/bash\n")
    lines.append("#SBATCH -J {}_{}\n".format(bt, file))
    lines.append("#SBATCH -n 6\n")
    lines.append("#SBATCH -N 1\n")
    lines.append("#SBATCH -p serial_requeue\n")
    lines.append("#SBATCH --mem 6000\n")
    lines.append("#SBATCH -t 00-10:00\n")
    lines.append("#SBATCH -o /n/home09/twagg/detecting-DCOs-in-LISA/simulation/slurm/logs/{}_{}_%A_%a.out\n".format(bt, file))
    lines.append("#SBATCH -e /n/home09/twagg/detecting-DCOs-in-LISA/simulation/slurm/logs/{}_{}_%A_%a.err\n".format(bt, file))
    lines.append("#SBATCH --mail-user thomas.wagg@cfa.harvard.edu\n")
    lines.append("#SBATCH --mail-type ALL\n")

    lines.append("module load Anaconda3/2020.11\n")
    lines.append("source activate lisa\n")

    main_line = 'python /n/home09/twagg/detecting-DCOs-in-LISA/simulation/src/simulate_DCO_detections.py -i /n/holystore01/LABS/berger_lab/Lab/fbroekgaarden/DATA/all_dco_legacy_CEbug_fix/{1}/COMPASOutputReduced.h5 -o /n/home09/twagg/detecting-DCOs-in-LISA/simulation/output/{0}_{1}_'.format(bt, file) + '"${SLURM_ARRAY_TASK_ID}".h5 -n 50 -t ' + '{}\n'.format(bt)
    if file == "optimistic":
        main_line = main_line.replace("\n", " --opt-flag\n")
    if simple_mw:
        main_line = 'python /n/home09/twagg/detecting-DCOs-in-LISA/simulation/src/simulate_DCO_detections.py -i /n/holystore01/LABS/berger_lab/Lab/fbroekgaarden/DATA/all_dco_legacy_CEbug_fix/{1}/COMPASOutputReduced.h5 -o /n/home09/twagg/detecting-DCOs-in-LISA/simulation/output/simple_mw_{0}_{1}_'.format(bt, file) + '"${SLURM_ARRAY_TASK_ID}".h5 -n 50 --simple-mw -t ' + '{}\n'.format(bt)
        script_file = script_file.replace(".sh", "_simplemw.sh")
        print(script_file)
    lines.append(main_line)

    
    with open(script_file, "w") as f:
        f.writelines(lines)


for bt in ["BHBH", "BHNS", "NSNS"]:
    create_file(bt, variations[0]["file"], simple_mw=True)
    for v in variations:
        create_file(bt, v["file"])
