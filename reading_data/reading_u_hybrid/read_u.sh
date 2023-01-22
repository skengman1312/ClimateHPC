#!/bin/bash
#PBS -N final_sub
#PBS -l select=2:ncpus=3:mpiprocs=5:mem=4gb 
#PBS -l walltime=0:09:00
#PBS -j oe
# set the excution on the short queue
#PBS -q short_cpuQ
cd ${PBS_O_WORKDIR}
module load Intel_parallel_studio_xe2019u1
module load mpich-3.2.1--intel-xe-2019u1
export I_MPI_DEBUG=5
export I_MPI_PIN_DOMAIN=omp
export I_MPI_PIN_ORDER=compact
export OMP_NUM_THREADS=5
export KMP_AFFINITY=verbose,compact

mpicc -std=c99 -g -Wall -fopenmp -I /apps/netCDF4.7.0--gcc-9.1.0/include -L /apps/netCDF4.7.0--gcc-9.1.0/lib -lnetcdf -o reading_u.out reading_u.c -lm -mcmodel medium
mpirun -n 5 $(pwd)/reading_u.out
date

ALLPROC=`cat $PBS_NODEFILE | wc -l`
echo The total number of nodes is $ALLPROC
export THIS_HOST=$(hostname)
echo Each place corresponds to a single $OMP_PLACES on the target machine. 
echo The host $THIS_HOST
