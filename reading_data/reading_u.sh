#!/bin/bash
#PBS -l select=5:ncpus=4:mem=4gb 

# set max execution time
#PBS -l walltime=0:10:00

# set the excution on the short queue
#PBS -q short_cpuQ
cd ${PBS_O_WORKDIR}
mpirun.actual -n 4 $(pwd)/reading_u.out
