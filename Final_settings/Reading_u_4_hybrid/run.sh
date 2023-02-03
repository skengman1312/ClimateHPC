#!/bin/bash
#PBS -l select=4:ncpus=4:mem=2gb 
#PBS -N run_reading_u_threaded
#PBS -j oe
#PBS -l walltime=0:20:00
#PBS -q short_cpuQ
cd ${PBS_O_WORKDIR} 
module load gcc91
module load mpich-3.2.1--gcc-9.1.0
# module load hwloc-2.0.4
module load openmpi-4.0.4
module load netcdf-4.7.0--gcc-9.1.0
module load hdf5-1.10.5--gcc-9.1.0
module list
export OMP_PLACES=threads
mpicc -g -Wall -I /apps/netCDF4.7.0--gcc-9.1.0/include -L /apps/netCDF4.7.0--gcc-9.1.0/lib -lnetcdf -o reading_u_hybrid.out reading_u_hybrid.c -lm -ldl -lz -lcurl -std=gnu99 -fopenmp
mpiexec -n 4 --report-bindings --map-by node:pe=4 --bind-to core  $(pwd)/reading_u_hybrid.out
# mpiexec -n 5  $(pwd)/reading_u_hybrid.out
date
ALLPROC=`cat $PBS_NODEFILE | wc -l`
echo The total number of nodes is $ALLPROC
export THIS_HOST=$(hostname)
## --map-by hwthread --bind-to core -genv OMP_PLACES=threads 
## --map-by node:pe=2