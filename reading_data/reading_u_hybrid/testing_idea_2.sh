#!/bin/bash
#PBS -l select=2:ncpus=2:mem=2gb 
#PBS -N testing_threading_mapping
#PBS -j oe
#PBS -l walltime=0:05:00
#PBS -q short_cpuQ
cd ${PBS_O_WORKDIR} 
module load gcc91
module load mpich-3.2.1--gcc-9.1.0
module load hwloc-2.0.4
module load openmpi-4.0.4
module load netcdf-4.7.0--gcc-9.1.0
module load hdf5-1.10.5--gcc-9.1.0
module list
export I_MPI_DEBUG=5
export MPI_MCA_mca_base_component_show_load_errors=0
export PMIX_MCA_mca_base_component_show_load_errors=0
ALLPROC=`cat $PBS_NODEFILE | wc -l`
echo The total number of nodes is $ALLPROC
export THIS_HOST=$(hostname)
export OMP_PLACES=threads
mpicc -g -Wall -I /apps/netCDF4.7.0--gcc-9.1.0/include -L /apps/netCDF4.7.0--gcc-9.1.0/lib -lnetcdf -o testing_idea.out testing_idea.c -lm -ldl -lz -lcurl -std=gnu99 -fopenmp
mpiexec -n 2 --report-bindings --map-by node:pe=2 --bind-to core  $(pwd)/testing_idea.out
date
## --map-by hwthread --bind-to core -genv OMP_PLACES=threads 
## --map-by node:pe=2