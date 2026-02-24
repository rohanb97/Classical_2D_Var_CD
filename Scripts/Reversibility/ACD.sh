#!/bin/bash -l

#$ -P f-dmrg
#$ -l h_rt=24:00:00
#$ -o scc_output/$JOB_NAME_$JOB_ID.txt
#$ -e scc_error/$JOB_NAME_$JOB_ID.txt

module load julia

julia scc_run/Energy_Fluc_bare.jl
