#!/bin/bash
#PBS -l select=20:ncpus=5:mem=1gb


# set max execution time
#PBS -l walltime=0:15:00


# set the excution on the short queue
#PBS -q short_cpuQ

cd ${PBS_O_WORKDIR}
mpicc -std=c99 -g -Wall -I /apps/netCDF4.7.0--gcc-9.1.0/include -L /apps/netCDF4.7.0--gcc-9.1.0/lib -lnetcdf -o reading_ssh_f.out reading_ssh_f.c -lm
module load mpich-3.2
mpirun.actual -n 20 $(pwd)/reading_ssh_f.out
