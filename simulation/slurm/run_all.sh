rm scripts/*.sh
python generate_scripts.py

for type in BHBH BHNS NSNS
do
    for variation in fiducial massTransferEfficiencyFixed_0_25 massTransferEfficiencyFixed_0_5 massTransferEfficiencyFixed_0_75 unstableCaseBB unstableCaseBB_opt alpha0_1 alpha0_5 alpha2_0 alpha10 optimistic rapid maxNSmass2_0 maxNSmass3_0 noPISN ccSNkick_100km_s ccSNkick_30km_s noBHkick wolf_rayet_multiplier_0_1 wolf_rayet_multiplier_5
    do
        sbatch --array=1-50 scripts/${type}_${variation}.sh
        sleep 1
    done
done