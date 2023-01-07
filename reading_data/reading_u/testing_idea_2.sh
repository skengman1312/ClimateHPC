#!/bin/bash
#PBS -l select=5:ncpus=2:ompthreads=2:mem=4gb 

# set max execution time
#PBS -l walltime=0:05:00

# set the excution on the short queue
#PBS -q short_cpuQ
cd ${PBS_O_WORKDIR}
mpirun -n 2 $(pwd)/testing_idea_2.out
