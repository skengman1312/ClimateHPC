#!/bin/bash
#PBS -l select=1:ncpus=4:mem=2gb

# set max execution time
#PBS -l walltime=0:05:00

# set the excution on the short queue
#PBS -q short_cpuQ
cd ${PBS_O_WORKDIR}

echo  pdw = $(pwd)/reading_data
module load mpich-3.2
mpirun.actual -n 4 $(pwd)/reading_data.out
