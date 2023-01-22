#!/bin/bash
#PBS -N hyb_4_3
#PBS -l select=5:ncpus=5:mpiprocs=5:ompthreads=5:mem=4gb 
#PBS -l walltime=0:10:00
#PBS -j oe
# set the excution on the short queue
#PBS -q short_cpuQ
cd ${PBS_O_WORKDIR} 
export I_MPI_DEBUG=5
ALLPROC=`cat $PBS_NODEFILE | wc -l`
echo The total number of nodes is $ALLPROC
echo MPIdebug is set to $I_MPI_DEBUG
export THIS_HOST=$(hostname)
export OMP_PLACES=threads
echo Each place corresponds to a single $OMP_PLACES on the target machine. 
echo The host $THIS_HOST
date
mpicc -std=c99 -g -Wall -fopenmp -I /apps/netCDF4.7.0--gcc-9.1.0/include -L /apps/netCDF4.7.0--gcc-9.1.0/lib -lnetcdf -o reading_u.out reading_u.c -lm 
mpirun -n 5 -genv I_MPI_DEBUG=5 --map-by hwthread --bind-to core $(pwd)/reading_u.out