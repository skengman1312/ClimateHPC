#!/bin/bash
#PBS -l select=5:ncpus=8:mem=2gb

# set max execution time
#PBS -l walltime=0:05:00

# set the excution on the short queue
#PBS -q short_cpuQ
cd ${PBS_O_WORKDIR}

module load mpich-3.2
mpirun.actual -n 8 /home/mostafa.haggag/programs/hello_world/mpi_peitro
