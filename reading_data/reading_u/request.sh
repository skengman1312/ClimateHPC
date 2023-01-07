#!/bin/bash
#PBS -l select=1:ncpus=4:mem=4gb 

# set max execution time
#PBS -l walltime=0:30:00 -I

# set the excution on the short queue
#PBS -q short_cpuQ
# cd ${PBS_O_WORKDIR}
