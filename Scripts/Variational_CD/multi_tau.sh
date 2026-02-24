#!/bin/bash -l

#$ -P f-dmrg
#$ -l h_rt=00:30:00
#$ -N E_Var_2D
#$ -o scc_output/$JOB_NAME_$JOB_ID.txt
#$ -e scc_error/$JOB_NAME_$JOB_ID.txt

#$ -v TAUPOINTS=10
#$ -v SAMPLEPOINTS=10

qsub -N ensemble -v MAIN_JOB_ID=$JOB_ID -v SAMPLE_POINTS=$SAMPLEPOINTS ensemble.sh

for order in `seq 0 1 1`; do
	qsub -N Var_Param_$order -hold_jid ensemble -v ORDER=$order -v MAIN_JOB_ID=$JOB_ID Var_Param.sh
	qsub -N ACD_$order -hold_jid Var_Param_$order -v ORDER=$order -v MAIN_JOB_ID=$JOB_ID -v POINTS=$TAUPOINTS -t 1-$TAUPOINTS ACD.sh
	qsub -N Clean_$order -hold_jid ACD_$order -v ORDER=$order -v MAIN_JOB_ID=$JOB_ID -v POINTS=$TAUPOINTS Combine.sh
done 	
