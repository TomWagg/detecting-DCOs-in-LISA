import sys
sys.path.append("../src/")
from variations import variations


def create_file(dco_type, file, simple_mw=False, extend_mission=False):
    extra = "_simplemw" if simple_mw else "_10yr" if extend_mission else ""
    script_file = "scripts/{}_{}{}.sh".format(dco_type, file, extra)

    lines = []
    lines.append("#!/bin/bash\n")
    lines.append("#SBATCH -J {}_{}{}\n".format(dco_type, file, extra))
    lines.append("#SBATCH -n 6\n")
    lines.append("#SBATCH -N 1\n")
    lines.append("#SBATCH -p serial_requeue\n")
    lines.append("#SBATCH --mem 6000\n")
    lines.append("#SBATCH -t 00-10:00\n")
    lines.append("#SBATCH -o /n/home09/twagg/detecting-DCOs-in-LISA/simulation/slurm/logs/{}_{}{}_%A_%a.out\n".format(dco_type, file, extra))
    lines.append("#SBATCH -e /n/home09/twagg/detecting-DCOs-in-LISA/simulation/slurm/logs/{}_{}{}_%A_%a.err\n".format(dco_type, file, extra))
    lines.append("#SBATCH --mail-user thomas.wagg@cfa.harvard.edu\n")
    lines.append("#SBATCH --mail-type ALL\n")

    lines.append("module load Anaconda3/2020.11\n")
    lines.append("source activate lisa\n")

    main_line = 'python /n/home09/twagg/detecting-DCOs-in-LISA/simulation/src/simulate_DCO_detections.py -i /n/holystore01/LABS/berger_lab/Lab/fbroekgaarden/DATA/all_dco_legacy_CEbug_fix/{1}/COMPASOutputCombined.h5 -o /n/home09/twagg/detecting-DCOs-in-LISA/simulation/output/{0}_{1}_'.format(dco_type, file) + '"${SLURM_ARRAY_TASK_ID}".h5 -n 50 -t ' + '{}\n'.format(dco_type)

    if simple_mw:
        main_line = 'python /n/home09/twagg/detecting-DCOs-in-LISA/simulation/src/simulate_DCO_detections.py -i /n/holystore01/LABS/berger_lab/Lab/fbroekgaarden/DATA/all_dco_legacy_CEbug_fix/{1}/COMPASOutputCombined.h5 -o /n/home09/twagg/detecting-DCOs-in-LISA/simulation/output/simple_mw_{0}_{1}_'.format(dco_type, file) + '"${SLURM_ARRAY_TASK_ID}".h5 -n 50 --simple-mw -t ' + '{}\n'.format(dco_type)
    elif extend_mission:
        main_line = 'python /n/home09/twagg/detecting-DCOs-in-LISA/simulation/src/simulate_DCO_detections.py -i /n/holystore01/LABS/berger_lab/Lab/fbroekgaarden/DATA/all_dco_legacy_CEbug_fix/{1}/COMPASOutputCombined.h5 -o /n/home09/twagg/detecting-DCOs-in-LISA/simulation/output/{0}_{1}_10yr_'.format(dco_type, file) + '"${SLURM_ARRAY_TASK_ID}".h5 -n 50 --extended-mission -t ' + '{}\n'.format(dco_type)

    if file == "optimistic":
        main_line = main_line.replace("\n", " --opt-flag\n").replace("all_dco_legacy_CEbug_fix/optimistic", "all_dco_legacy_CEbug_fix/fiducial")
    if file == "unstableCaseBB_opt":
        main_line = main_line.replace("\n", " --opt-flag\n").replace("all_dco_legacy_CEbug_fix/unstableCaseBB_opt", "all_dco_legacy_CEbug_fix/unstableCaseBB")
    lines.append(main_line)


    with open(script_file, "w") as f:
        f.writelines(lines)


for dco_type in ["BHBH", "BHNS", "NSNS"]:
    create_file(dco_type, "fiducial", simple_mw=True)
    create_file(dco_type, "unstableCaseBB_opt", simple_mw=True)

    for extend_mission in [False, True]:
        for v in variations:
            create_file(dco_type, v["file"], extend_mission=extend_mission)
