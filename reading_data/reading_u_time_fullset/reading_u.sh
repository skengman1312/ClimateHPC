#!/bin/bash
#PBS -l select=1:ncpus=4:ompthreads=1:mem=4gb 
#PBS -j eo
# set max execution time
#PBS -l walltime=0:30:00

# set the excution on the short queue
#PBS -q short_cpuQ
cd ${PBS_O_WORKDIR}
mpicc -std=c99 -g -Wall -fopenmp -I /apps/netCDF4.7.0--gcc-9.1.0/include -L /apps/netCDF4.7.0--gcc-9.1.0/lib -lnetcdf -o reading_u.out reading_u.c -lm 
export I_MPI_DEBUG=5
mpirun -np 4  $(pwd)/reading_u.out
