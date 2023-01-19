#!/bin/bash
#PBS -N hyb_4_3
#PBS -l select=1:ncpus=5:mpiprocs=1:ompthreads=5:mem=4gb 
#PBS -l walltime=0:05:00
#PBS -j oe
# set the excution on the short queue
#PBS -q short_cpuQ
module load hwloc-2.0.4
module load Intel_parallel_studio_xe2019u1
module load mpich-3.2.1--intel-xe-2019u1
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
mpicc -std=c99 -g -Wall -fopenmp -I /apps/netCDF4.7.0--gcc-9.1.0/include -L /apps/netCDF4.7.0--gcc-9.1.0/lib -lnetcdf -o file.out file.c -lm 
mpirun -n 5 -genv I_MPI_DEBUG=5 --map-by hwthread --bind-to core $(pwd)/file.out