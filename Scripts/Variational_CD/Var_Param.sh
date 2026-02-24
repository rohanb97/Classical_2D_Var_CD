#!/bin/bash -l

#$ -P f-dmrg
#$ -l h_rt=02:00:00
#$ -o scc_output/$JOB_NAME_$JOB_ID.txt
#$ -e scc_error/$JOB_NAME_$JOB_ID.txt

module load julia

julia scc_run/Compute_Var_Param.jl
