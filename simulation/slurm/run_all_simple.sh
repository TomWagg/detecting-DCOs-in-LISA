rm scripts/*.sh
python generate_scripts.py

for type in BHBH BHNS NSNS
do
    for variation in fiducial unstableCaseBB_opt
    do
        sbatch --array=1-50 scripts/${type}_${variation}_simplemw.sh
        sleep 1
    done
done