#!/bin/bash
#PBS -l select=2:ncpus=1:mem=2gb 

# set max execution time
#PBS -l walltime=0:10:00

# set the excution on the short queue
#PBS -q short_cpuQ
cd ${PBS_O_WORKDIR}
mpirun.actual -n 4 $(pwd)/reading_u
