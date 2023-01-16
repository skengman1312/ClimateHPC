#!/bin/bash
#PBS -l select=3:ncpus=4:mem=4gb 

# set max execution time
#PBS -l walltime=0:10:00

# set the excution on the short queue
#PBS -q short_cpuQ
cd ${PBS_O_WORKDIR}
mpicc -std=c99 -g -Wall -I /apps/netCDF4.7.0--gcc-9.1.0/include -L /apps/netCDF4.7.0--gcc-9.1.0/lib -lnetcdf -o reading_u.out reading_u.c -lm 

mpirun -n 9 $(pwd)/reading_u.out
