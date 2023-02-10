#!/bin/bash
#PBS -N run_reading_u_simple
#PBS -l select=4:ncpus=4:mem=8gb 
#PBS -j oe
# set max execution time
#PBS -l walltime=0:30:00

# set the excution on the short queue
#PBS -q short_cpuQ


cd ${PBS_O_WORKDIR}
module load gcc91
module load mpich-3.2.1--gcc-9.1.0
module load netcdf-4.7.0--gcc-9.1.0
module load hdf5-1.10.5--gcc-9.1.0
module list
echo "Starting the serial version" 
mpicc -std=c99 -g -Wall -fopenmp -I /apps/netCDF4.7.0--gcc-9.1.0/include -L /apps/netCDF4.7.0--gcc-9.1.0/lib -lnetcdf -o serial.out serial.c -lm 
mpirun -np 1 -prepend-rank $(pwd)/serial.out
echo "Ending the serial version" 

echo "Starting the simple version" 
mpicc -std=c99 -g -Wall -fopenmp -I /apps/netCDF4.7.0--gcc-9.1.0/include -L /apps/netCDF4.7.0--gcc-9.1.0/lib -lnetcdf -o reading_u.out reading_u_simple.c -lm 
mpirun -np 16 -prepend-rank $(pwd)/reading_u.out

echo "Ending the simple version"

echo "Starting the comm splitting version" 
mpicc -g -Wall -I /apps/netCDF4.7.0--gcc-9.1.0/include -L /apps/netCDF4.7.0--gcc-9.1.0/lib -lnetcdf -o reading_u_hybrid.out reading_u_spill_comm_v2.c -lm -ldl -lz -lcurl -std=gnu99 -fopenmp
mpirun -np 16 -prepend-rank  $(pwd)/reading_u_hybrid.out
echo "Ending the comm splitting version"
export OMP_PLACES=threads
module load openmpi-4.0.4

echo "Starting the comm splitting version" 
mpicc -g -Wall -I /apps/netCDF4.7.0--gcc-9.1.0/include -L /apps/netCDF4.7.0--gcc-9.1.0/lib -lnetcdf -o reading_u_hybrid_v1.out reading_u_hybrid_v1.c -lm -ldl -lz -lcurl -std=gnu99 -fopenmp
mpiexec -n 4 --report-bindings --map-by node:pe=4 --bind-to core  $(pwd)/reading_u_hybrid_v1.out
echo "Ending the comm splitting version"

echo "Starting the comm splitting version" 
mpicc -g -Wall -I /apps/netCDF4.7.0--gcc-9.1.0/include -L /apps/netCDF4.7.0--gcc-9.1.0/lib -lnetcdf -o reading_u_hybrid_v2.out reading_u_hybrid_v2.c -lm -ldl -lz -lcurl -std=gnu99 -fopenmp
mpiexec -n 4 --report-bindings --map-by node:pe=4 --bind-to core  $(pwd)/reading_u_hybrid_v2.out

echo "Ending the comm splitting version"
date

ALLPROC=`cat $PBS_NODEFILE | wc -l`
echo The total number of nodes is $ALLPROC
export THIS_HOST=$(hostname)
echo Each place corresponds to a single $OMP_PLACES on the target machine. 
echo The host $THIS_HOST