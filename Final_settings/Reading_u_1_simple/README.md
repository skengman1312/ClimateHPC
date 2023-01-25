# Toplogy
* The main idea in here is to make eeach process read serially a part of the file. 
* You donot use any threading in here.
* You only use processses.
* The idea is not to use many processes and divide the data. 

mpicc -std=c99 -g -Wall -fopenmp -I /apps/netCDF4.7.0--gcc-9.1.0/include -L /apps/netCDF4.7.0--gcc-9.1.0/lib -lnetcdf -o reading_u.out reading_u_simple.c -lm 
