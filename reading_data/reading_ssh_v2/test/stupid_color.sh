#!/bin/bash
#PBS -l select=1:ncpus=12:mem=2gb

# set max execution time
#PBS -l walltime=0:10:00


# set the excution on the short queue
#PBS -q short_cpuQ

cd ${PBS_O_WORKDIR}
mpicc -std=c99 -g -Wall -I /apps/netCDF4.7.0--gcc-9.1.0/include -L /apps/netCDF4.7.0--gcc-9.1.0/lib -lnetcdf -o stupid_color.out stupid_color.c -lm
module load mpich-3.2
mpirun.actual -n 12 $(pwd)/stupid_color.out
