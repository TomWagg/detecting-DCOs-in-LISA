rm scripts/*.sh
python generate_scripts.py

for type in BHBH BHNS NSNS
do
    for variation in fiducial optimistic alpha0_5 alpha2_0 unstableCaseBB massTransferEfficiencyFixed_0_25 massTransferEfficiencyFixed_0_5 massTransferEfficiencyFixed_0_75 ccSNkick_100km_s ccSNkick_30km_s noBHkick rapid maxNSmass2_0 maxNSmass3_0 noPISN
    do
        sbatch --array=1-50 scripts/${type}_${variation}.sh
        sleep 10
    done
done