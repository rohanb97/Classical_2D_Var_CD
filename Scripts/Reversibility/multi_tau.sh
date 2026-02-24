#!/bin/bash -l
#$ -P f-dmrg
#$ -l h_rt=00:30:00
#$ -N E_Var_2D
#$ -o scc_output/$JOB_NAME_$JOB_ID.txt
#$ -e scc_error/$JOB_NAME_$JOB_ID.txt
#$ -v SAMPLEPOINTS=20000

# Read tau and beta lists from config file
CONFIG_FILE="$HOME/classical-cd/bear-drive/Julia_Files/tau_beta_config.txt"

# Extract non-comment lines
TAULIST_STR=$(grep -v '^#' "$CONFIG_FILE" | sed -n '1p' | tr -d ' ')
BETALIST_STR=$(grep -v '^#' "$CONFIG_FILE" | sed -n '2p' | tr -d ' ')

# Convert comma-separated strings to bash arrays
IFS=',' read -ra TAULIST <<< "$TAULIST_STR"
IFS=',' read -ra BETALIST <<< "$BETALIST_STR"

TAUPOINTS=${#TAULIST[@]}
BETAPOINTS=${#BETALIST[@]}

echo "Read from config file:"
echo "Tau values: ${TAULIST[@]}"
echo "Beta values: ${BETALIST[@]}"

# Submit ensemble job
qsub -N ensemble -v MAIN_JOB_ID=$JOB_ID -v SAMPLE_POINTS=$SAMPLEPOINTS ensemble.sh

# Submit jobs for each tau-beta combination
task_id=1
for ((i=0; i<$TAUPOINTS; i++)); do
    for ((j=0; j<$BETAPOINTS; j++)); do
        TAU=${TAULIST[$i]}
        BETA=${BETALIST[$j]}
        qsub -N ACD_${i}_${j} -hold_jid ensemble \
             -v MAIN_JOB_ID=$JOB_ID \
             -v TAU_VALUE=$TAU -v BETA_VALUE=$BETA \
             -v TASK_ID=$task_id ACD.sh
        ((task_id++))
    done
done

# Submit cleanup job
qsub -N Clean -hold_jid "ACD_*" -v MAIN_JOB_ID=$JOB_ID Combine.sh
