#!/bin/bash
#SBATCH --mail-user=nguyenk@ufl.edu
#SBATCH --job-name=sofetch
#SBATCH --mail-type=FAIL
#SBATCH --output=/ufrc/ewang/nguyenk/nmd_reg/outlog_fetch/nmd_reg_%j.out
#SBATCH --error=/ufrc/ewang/nguyenk/nmd_reg/errorlog_fetch/nmd_reg_%j.out
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16 ## threads
#SBATCH --nodes=1
#SBATCH --mem=5gb
#SBATCH --array=1
-2999%20
#SBATCH --qos=ewang-b

module load python

PROJ=/ufrc/ewang/nguyenk/nmd_reg

DESIGN_FILE=${PROJ}/design_rnalist_submit1.csv
DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "${DESIGN}"

RNA=${ARRAY[0]}

DATA="${PROJ}/human.txt"

echo ${RNA}

python make_bed.py --rna ${RNA} --data ${DATA}

ID=$(<id_${RNA}.txt)

echo "Made BED Files for ${RNA}. ID is ${ID}. Directory is ${RNA}bed${ID}/"

BED_DIR="${PROJ}/${RNA}bed${ID}/*"

mkdir "${PROJ}/${RNA}made_beds${ID}/"
echo "Made directory for fetchbeds. Name is ${RNA}made_beds${ID}/."

for BED in ${BED_DIR}
do
	TEMPID=$(echo $BED| cut -d'_' -f 4)
	echo ${TEMP}
    python fetchFromBed.py --fetch ${BED} chroms/ "${PROJ}/${RNA}made_beds${ID}/${RNA}_bed_${TEMPID}_made.txt"

    echo "Completed fetch for ${BED#*/}."
done

echo "Starting startstop process."
python startstop.py --rna ${RNA} --id ${ID}
echo "Startstop process complete, output to ${CSV_DIR}ss_${RNA}-${ID}.csv."

echo "Starting setflags process."
python setflags.py --dir ${PROJ} --rna ${RNA} --id ${ID}
echo "Setflags process complete, output to ${CSV_DIR}flag_${RNA}-${ID}.csv."

CSV_DIR="${PROJ}/${RNA}-${ID}_analysis/"
mkdir ${CSV_DIR}

mv ss_${RNA}-${ID}.csv ${CSV_DIR}
mv flag_${RNA}-${ID}.csv ${CSV_DIR}

FINAL_DIR="${PROJ}/${RNA}-${ID}_all_data/"
mkdir ${FINAL_DIR}

mv id_${RNA}.txt ${FINAL_DIR}
mv ${RNA}bed${ID} ${FINAL_DIR}
mv ${RNA}bedinfo${ID} ${FINAL_DIR}
mv ${RNA}made_beds${ID} ${FINAL_DIR}
mv ${RNA}-${ID}_analysis ${FINAL_DIR}

mv ${FINAL_DIR} ${PROJ}/data_dir/

echo "Bash script complete. All information outputted into ${PROJ}/data_dir/${FINAL_DIR}."