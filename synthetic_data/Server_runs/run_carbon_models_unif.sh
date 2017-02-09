#!/bin/sh

# Directives
#PBS -N power_analysis_Carbon
#PBS -W group_list=yetistats
#PBS -l nodes=1:ppn=2:v2,walltime=24:00:00,mem=8000mb
#PBS -t 1-18
#PBS -M n 
#PBS -m abe
#PBS -V

# Set output and error directories
#PBS -o localhost:/vega/stats/users/mk3971/Carbon/server_runs/outputs/unif_
#PBS -e localhost:/vega/stats/users/mk3971/Carbon/server_runs/outputs/unif_

# Run rscript on each node
Rscript simulate_and_fit_uniform.R

# End of script
