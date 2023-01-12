#!/bin/bash
#PBS -l select=2:ncpus=4:ompthreads=8:mem=4gb 

# set max execution time
#PBS -l walltime=0:15:00

# set the excution on the short queue
#PBS -q short_cpuQ
cd ${PBS_O_WORKDIR}
mpicc -std=c99 -g -Wall -fopenmp -I /apps/netCDF4.7.0--gcc-9.1.0/include -L /apps/netCDF4.7.0--gcc-9.1.0/lib -lnetcdf -o reading_v.out reading_v.c -lm
mpirun.actual -n 8 $(pwd)/reading_v.out
