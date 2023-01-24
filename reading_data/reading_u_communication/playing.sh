#!/bin/bash
#PBS -N playing
#PBS -l select=1:ncpus=2:mem=10gb 
#PBS -l walltime=0:1:00
#PBS -j oe
# set the excution on the short queue
#PBS -q short_cpuQ
cd ${PBS_O_WORKDIR}

export I_MPI_DEBUG=5
mpicc -std=c99 -g -Wall -fopenmp -I /apps/netCDF4.7.0--gcc-9.1.0/include -L /apps/netCDF4.7.0--gcc-9.1.0/lib -lnetcdf -o playing.out playing.c -lm 
mpirun -np 2 ./playing.out
date

ALLPROC=`cat $PBS_NODEFILE | wc -l`
echo The total number of nodes is $ALLPROC
export THIS_HOST=$(hostname)
echo Each place corresponds to a single $OMP_PLACES on the target machine. 
echo The host $THIS_HOST
## -mcmodel medium