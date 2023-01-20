#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <netcdf.h>
#include <mpi.h>
#include <omp.h>
#include <sys/time.h>
/*
mpicc -std=c99 -g -Wall -fopenmp -I /apps/netCDF4.7.0--gcc-9.1.0/include -L /apps/netCDF4.7.0--gcc-9.1.0/lib -lnetcdf -o reading_u.out reading_u.c -lm 
*/
/*MACROS START*/
#define FILE_NAME "/shares/HPC4DataScience/FESOM2/unod.fesom.2010.nc"
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); return 2;} //just simple macro to handle error with NETCDF
#define HEIGHT_1_NAME "nz1"  // the variable is called nz in the file. I knew it using ncdump command
#define TIME "time"  // the variable is called time in the file. I knew it using ncdump command
#define UNOD "unod"  // the variable is called unod in the file. I knew it using ncdump command
#define N_NZ1 69            // the variable nz1 has 69 entries
#define N_TIME 12            // the variable time has 12 entries
#define NDIMS 3
#define GRID_POINTS 8852366
#define DEBUG 1
// for write
#define FILE_NAME2 "map_summarized_fuldataset.nc"
#define UNITS_speed "m_s"
#define UNITS "units"
#define UNITS_time "s"
#define NDIMS_wr 3
/*MACROS end*/

/*FUNCTIONS START*/
/*This is a function that measures time using system time val
We subtract both time instances and we convert final output in seconds*/
double time_diff(struct timeval *start, struct timeval *end);
void convert_time_hour_sec( double seconds, long int *h, long int *m, long int *s);
/*FUNCTIONS end*/


int main (int argc, char *argv[]){
    /* MPI  inizialization */
    MPI_Init(&argc, &argv);

    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);
    // Print off a hello world message with processor name
    printf("greetings:  %s, rank %d out of %d processes\n",
           processor_name, rank, size);

    /*NETCDF id*/
    int ncid;
    int unod_id;
    int retval;
    /*Netcdf id for write*/
    // int ncid2; // write file nmae
    // int time_dimid; // time dimension id 
    // int speed_dimid; // speed dimension id 
    // int var_speed_id;// variable speed id 
    size_t start[NDIMS], count[NDIMS];
    /*LOOOPING variables*/
    // int rec;
    int i;
    int k ;
    int levels_per_proc = ceil((double)N_NZ1 / size);/*if 2.3 is passed to ceil(), it will return 3*/
    float levels[levels_per_proc][GRID_POINTS];
    int sendcounts[size];// for scatterv
    int displs[size];// for scatterv
    // intializing the arrays 
    for (i = 0; i < size;i++){
        sendcounts[i]= levels_per_proc * GRID_POINTS;
        displs[i] = i * levels_per_proc * GRID_POINTS;
    }
    sendcounts[size-1] = (N_NZ1-(levels_per_proc*(size-1)))*GRID_POINTS;
    static float u_speed[N_NZ1][GRID_POINTS] = {{0}};
    if(0==rank){
        if ((retval = nc_open(FILE_NAME, NC_NOWRITE, &ncid)))
        ERR(retval);
        if ((retval = nc_inq_varid(ncid, UNOD, &unod_id)))
        ERR(retval);
        count[0] = 1;/*1 time*/
        count[1] = N_NZ1;/*1 level*/
        count[2] = GRID_POINTS;/*all gridpoints*/
        start[1] = 0;
        start[2] = 0;
    }
    /*START*/
    for (k = 0;k<N_TIME; k++){
        if(rank==0){
        // int sendcounts[size];// fill it all with 17 *GRID_POINTS except the last one 15*GRID_POINTS
            start[0] = k;
            if ((retval = nc_get_vara_float(ncid, unod_id, start, 
                        count, &u_speed[0][0])))
                        ERR(retval);
            // printf("for process 0 =%.6f\n", u_speed[0][8852365]);
            // printf("for process 0 =%.6f\n", u_speed[14][8852365]);
            // printf("for process 0 =%.6f\n", u_speed[28][8852365]);
            // printf("for process 0 =%.6f\n", u_speed[42][8852365]);
            // printf("for process 0 =%.6f\n", u_speed[56][8852365]);
            MPI_Scatterv(u_speed,sendcounts,displs,MPI_FLOAT,levels,GRID_POINTS*levels_per_proc,MPI_FLOAT,0,MPI_COMM_WORLD);
            // printf("%.6f from process 0\n", levels[0][8852365]);
        }
        else{
            MPI_Scatterv(u_speed,sendcounts,displs,MPI_FLOAT,levels,GRID_POINTS*levels_per_proc,MPI_FLOAT,0,MPI_COMM_WORLD);
            // printf("%.6f from process %d\n", levels[0][8852365],rank);
            }
    }

   MPI_Finalize();

   return 0;

}



double time_diff(struct timeval *start, struct timeval *end)
{ 
    // long tv_sec;                /* seconds */
    // long tv_usec;               /* microseconds */ 
    return ((end->tv_sec - start->tv_sec) +  1e-6 *(end->tv_usec - start->tv_usec));
}
void convert_time_hour_sec(double seconds,long int *h,long int *m,long int *s)
{ 	
    *h = ((long int) seconds/3600); 
	*m = ((long int) seconds -(3600 * ( *h )))/60;
    *s = ((long int) seconds -(3600 * ( *h ))-(( *m ) * 60));
}