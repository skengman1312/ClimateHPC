#!/bin/bash
#PBS -l select=2:ncpus=2:ompthreads=4:mem=4gb 

# set max execution time
#PBS -l walltime=0:30:00

# set the excution on the short queue
#PBS -q short_cpuQ
cd ${PBS_O_WORKDIR}
module load openmpi-4.0.4
export  OMP_DISPLAY_ENV="TRUE" 
mpicc -std=c99 -g -Wall -fopenmp -I /apps/netCDF4.7.0--gcc-9.1.0/include -L /apps/netCDF4.7.0--gcc-9.1.0/lib -lnetcdf -o reading_u.out reading_u.c -lm 

mpirun -n 4 $(pwd)/reading_u.out
