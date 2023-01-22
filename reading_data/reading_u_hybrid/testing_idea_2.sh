#!/bin/bash
#PBS -l select=1:ncpus=5:ompthreads=1:mem=2gb 
#PBS -N testing_threading_mapping
#PBS -j oe
#PBS -l walltime=0:05:00
#PBS -q short_cpuQ
cd ${PBS_O_WORKDIR}
export I_MPI_DEBUG=5
ALLPROC=`cat $PBS_NODEFILE | wc -l`
echo The total number of nodes is $ALLPROC
export THIS_HOST=$(hostname)
export OMP_PLACES=threads
mpicc -std=c99 -g -Wall -fopenmp -I /apps/netCDF4.7.0--gcc-9.1.0/include -L /apps/netCDF4.7.0--gcc-9.1.0/lib -lnetcdf -o testing_idea.out testing_idea.c -lm 
mpirun -n 1  --map-by hwthread --bind-to core -genv OMP_PLACES=threads  $(pwd)/testing_idea.out
date
## --map-by hwthread --bind-to core -genv OMP_PLACES=threads 
