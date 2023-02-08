#!/bin/bash
#PBS -N run_reading_u_simple
#PBS -l select=1:ncpus=16:ompthreads=1:mem=4gb 
#PBS -j oe
# set max execution time
#PBS -l walltime=0:20:00

# set the excution on the short queue
#PBS -q short_cpuQ


cd ${PBS_O_WORKDIR}
module load gcc91
module load mpich-3.2.1--gcc-9.1.0
module load netcdf-4.7.0--gcc-9.1.0
module load hdf5-1.10.5--gcc-9.1.0
module list

export I_MPI_DEBUG=5

mpicc -std=c99 -g -Wall -fopenmp -I /apps/netCDF4.7.0--gcc-9.1.0/include -L /apps/netCDF4.7.0--gcc-9.1.0/lib -lnetcdf -o reading_u.out reading_u_simple.c -lm 

mpirun -np 16 -prepend-rank $(pwd)/reading_u.out

date

ALLPROC=`cat $PBS_NODEFILE | wc -l`
echo The total number of nodes is $ALLPROC
export THIS_HOST=$(hostname)
echo Each place corresponds to a single $OMP_PLACES on the target machine. 
echo The host $THIS_HOST