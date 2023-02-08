#!/bin/bash
#PBS -l select=1:ncpus=1:mem=2gb

# set max execution time
#PBS -l walltime=0:0:05


# set the excution on the short queue
#PBS -q short_cpuQ

cd ${PBS_O_WORKDIR}
mpicc -std=c99 -g -Wall -I /apps/netCDF4.7.0--gcc-9.1.0/include -L /apps/netCDF4.7.0--gcc-9.1.0/lib -lnetcdf -o prova.out prova.c -lm
module load mpich-3.2
mpirun.actual -n 1 $(pwd)/prova.out
