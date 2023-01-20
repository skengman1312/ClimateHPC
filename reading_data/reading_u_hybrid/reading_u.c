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


/*MACROS end*/

/*FUNCTIONS START*/
/*This is a function that measures time using system time val
We subtract both time instances and we convert final output in seconds*/
double time_diff(struct timeval *start, struct timeval *end);
void convert_time_hour_sec( double seconds, long int *h, long int *m, long int *s);
void prepare_file(int * file_id,int * speed_id);
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
    /* VARIABLES DEFINE START*/
    /*
    
    * The is idea is that you loop over time 
    * You read big data
    * You scatter on all of the rest each has some amount of levels
    * you recieved  18 level * 8 million
    * you thread over 2 functions 
    * you gather everything 
    
    */

    /*NETCDF id*/
    int ncid;
    int unod_id;
    int retval;
    /*Netcdf id for write*/
    int ncid2; // write file nmae
    int time_dimid; // time dimension id 
    int speed_dimid; // speed dimension id 
    int var_speed_id;// variable speed id 
    size_t start[NDIMS], count[NDIMS];
    /*LOOOPING variables*/
    int rec;
    int i;
    int k ;
    static float x[11][GRID_POINTS];
    /*START*/
    if(0==rank){
        prepare_file(&ncid,&unod_id);
        int levels_per_proc = ceil((double)N_NZ1 / size);/*if 2.3 is passed to ceil(), it will return 3*/
        // int limit = (rank + 1) * levels_per_proc;

        count[0] = 1;/*1 time*/
        count[1] = N_NZ1;/*1 level*/
        count[2] = GRID_POINTS;/*all gridpoints*/
        // start[0] = k;
        start[1] = 0;
        start[2] = 0;
        int sendcounts[size];// fill it all with 17 *GRID_POINTS except the last one 15*GRID_POINTS
        int displs[size];// fill it with zeros
        float **u_speed = malloc(sizeof(float *) * N_NZ1);
        for (int i = 0; i < N_NZ1; i++){
            u_speed[i] = (float *)calloc(GRID_POINTS ,sizeof(float));
        }
        for (k = 0;k<N_TIME; k++){
                start[0] = k;
            if ((retval = nc_get_vara_float(ncid, unod_id, start, 
                        count, &u_speed[0][0])))
                        ERR(retval);
            MPI_Scatterv(u_speed,sendcounts,displs,MPI_FLOAT,x,sendcounts,MPI_FLOAT,0,MPI_COMM_WORLD);

        }
    }
    else{

    }
    
}
void prepare_file(int * file_id,int * speed_id){
    int retval;
    if ((retval = nc_open(FILE_NAME, NC_NOWRITE, &file_id)))
        ERR(retval);
    if ((retval = nc_inq_varid(file_id, UNOD, &speed_id)))
        ERR(retval);
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