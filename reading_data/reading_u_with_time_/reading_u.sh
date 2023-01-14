#!/bin/bash
#PBS -l select=2:ncpus=4:ompthreads=8:mem=4gb 

# set max execution time
#PBS -l walltime=0:30:00

# set the excution on the short queue
#PBS -q short_cpuQ
cd ${PBS_O_WORKDIR}
mpirun -n 8 $(pwd)/reading_u.out
