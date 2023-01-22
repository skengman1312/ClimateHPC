#!/bin/bash
#PBS -l select=1:ncpus=5:mpiprocs=1:mem=2gb 
#PBS -N threads_baseline
#PBS -j oe
#PBS -l walltime=0:05:00
#PBS -q short_cpuQ
cd ${PBS_O_WORKDIR}
ALLPROC=`cat $PBS_NODEFILE | wc -l`
echo The total number of nodes is $ALLPROC
module load Intel_parallel_studio_xe2019u1
module load mpich-3.2.1--intel-xe-2019u1
export I_MPI_DEBUG=5
export THIS_HOST=$(hostname)
export I_MPI_DEBUG=5
export I_MPI_PIN_DOMAIN=omp
export I_MPI_PIN_ORDER=compact
export OMP_NUM_THREADS=1
export KMP_AFFINITY=verbose,compact

mpicc -std=c99 -g -Wall -fopenmp -I /apps/netCDF4.7.0--gcc-9.1.0/include -L /apps/netCDF4.7.0--gcc-9.1.0/lib -lnetcdf -o testing_idea.out testing_idea.c -lm 
mpirun -prepend-rank   $(pwd)/testing_idea.out
date
## --map-by hwthread --bind-to core -genv OMP_PLACES=threads 
## pid Show process id for each debug message.
## tid	Show thread id for each debug message for multithreaded library.