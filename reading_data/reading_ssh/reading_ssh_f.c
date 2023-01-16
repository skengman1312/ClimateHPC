//
// Created by Pietro on 16/01/2023.
//
#include <stdio.h>
#include <string.h>
#include <netcdf.h>
#include<mpi.h>
/*
mpicc -std=c99 -g -Wall -I /apps/netCDF4.7.0--gcc-9.1.0/include -L /apps/netCDF4.7.0--gcc-9.1.0/lib -lnetcdf -o reading_ssh.out reading_ssh.c -lm
*/
/*MACROS START*/
#define FILE_NAME "/shares/HPC4DataScience/ssh.fesom.2010.nc"
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); return 2;} //just simple macro to handle error with NETCDF
#define HEIGHT_1_NAME "nz1"  // the variable is called nz in the file. I knew it using ncdump command
#define TIME "time"  // the variable is called time in the file. I knew it using ncdump command
#define UNOD "vnod"  // the variable is called vnod in the file. I knew it using ncdump command
#define SSH "ssh"  // the variable is called vnod in the file. I knew it using ncdump command
#define N_NZ1 69            // the variable nz1 has 69 entries
#define N_TIME 365            // the variable time has 12 entries
#define NDIMS 2
#define GRID_POINTS 8852366

// for write
#define FILE_NAME2 "map_summarized.nc"
#define UNITS_ssh "m"
#define UNITS "units"
#define UNITS_time "s"

/*MACROS end*/

float * average(float local_ssh[][GRID_POINTS], int timeframe);

int main () {
    // Initialize the MPI environment. The two arguments to MPI Init are not
    // currently used by MPI implementations, but are there in case future
    // implementations might need the arguments.
    MPI_Init(NULL, NULL);

    // /* VARIABLES DEFINE START*/
    // get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // get the rank of processes
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

//    // size across each dimension.
//    int dim[2] = {30, GRID_POINTS};
//    // dimensions of the scattered data received by each process
//    int local_dim[2] = {dim[0] / world_size, dim[1]};
//
//    // receive buffer
//    float rec[local_dim[0]][local_dim[1]];
//    float loc_avg[local_dim[1]];
//    float avg[local_dim[1]];
//    int sendcnt = local_dim[0] * local_dim[1]; /* how many items are sent to each process */
//    int recvcnt = local_dim[0] * local_dim[1];

    /*NETCDF id*/
    int ncid;
    int ssh_varid;
    int time_varid;
    int ssh_id;
    int retval;

    /*Netcdf id for write*/
    int ncid2;
    int time_new_id;
    int gp_new_id;
    int time_var_new_id;
    int var_new_id;



    /* Program variables to hold the data we will read.
    We will need enough space to hold all the timesteps of data. */
    static float ssh[N_TIME][GRID_POINTS];

    /* The start and count arrays will tell the netCDF library where to read our data. */
    size_t start[NDIMS], count[NDIMS];
    if (world_rank == 0) {
        /* Open the file. */
        if ((retval = nc_open(FILE_NAME, NC_NOWRITE, &ncid))) ERR(retval);

        /* Get the varids of the time and nz1 coordinate variables. */
        if ((retval = nc_inq_varid(ncid, SSH, &ssh_varid))) ERR(retval);
        if ((retval = nc_inq_varid(ncid, TIME, &time_varid))) ERR(retval);


        /* Get the varid of the ssh netCDF variable. */
        if ((retval = nc_inq_varid(ncid, SSH, &ssh_id))) ERR(retval);
        /* Read the data. Since we know the contents of the file we know that the
        data arrays in this program are the correct size to hold one timestep.*/
        count[0] = 1;
//    count[1] = N_NZ1;
        count[1] = GRID_POINTS;
        start[0] = 0;
        start[1] = 0;
        /* end of setup of NetCDF reading */

        /*LOOOPING variables*/

        /* sum matrix */
//    static float sum_u_speed[N_NZ1][GRID_POINTS] = {{0}};
        for (int i = 0; i < 30+30; i++) {
            start[0] = i;
            if ((retval = nc_get_vara_float(ncid, ssh_id, start,
                                            count, &ssh[i][0]))) ERR(retval);

        }
        printf("%lf\n", ssh[0+30][0]);
        printf("%lf\n", ssh[29+30][8852365]);
        /*CLOSING FILE*/
        if ((retval = nc_close(ncid))) ERR(retval);
        printf("*** SUCCESS :) reading example %s\n", FILE_NAME);
        // End of I

    }



    average(ssh+30, 30);

    MPI_Finalize();
    return 0;
}

float * average(float local_ssh[][GRID_POINTS], int timeframe ){
    // /* VARIABLES DEFINE START*/
    // get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // get the rank of processes
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // size across each dimension.
    int dim[2] = {timeframe, GRID_POINTS};
    // dimensions of the scattered data received by each process
    int local_dim[2] = {dim[0] / world_size, dim[1]};

    // receive buffer
    float rec[local_dim[0]][local_dim[1]];
    float loc_avg[local_dim[1]];
    float *avg = malloc(sizeof(local_dim[1])
    float avg[local_dim[1]];
    int sendcnt = local_dim[0] * local_dim[1]; /* how many items are sent to each process */
    int recvcnt = local_dim[0] * local_dim[1];

    printf("local ssh first element: %lf\n", local_ssh[0][0]);
    // MPI call and functional part
    MPI_Scatter(local_ssh, sendcnt, MPI_FLOAT,
                rec, recvcnt, MPI_FLOAT, 0, MPI_COMM_WORLD);
    printf("My rank is %i\n", world_rank);
    // loading the sum in a vector of size N-points
    for (int i = 0; i < local_dim[1]; i++){
        for (int j = 0; j < local_dim[0]; j++) {
            loc_avg[i] += rec[j][i];
        }
    }

    // averaging the sums for every element in the vector
    for (int i = 0; i < local_dim[1]; ++i) {
        loc_avg[i] = loc_avg[i] / local_dim[0];
    }

    printf("My rank is %i and I received from  %g to %g\n", world_rank, rec[0][0],rec[local_dim[0]-1][local_dim[1]-1]);
    printf("My loc_avg array is %g, %g, %g\n", loc_avg[0], loc_avg[1], loc_avg[2]);

    // reduce with sum operator
    MPI_Reduce(loc_avg, avg, local_dim[1], MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

    // averaging the results of the several threads to get as single average vector
    for (int i = 0; i < local_dim[1]; ++i) {
        avg[i] = avg[i]/world_size;
    }
    if (world_rank == 0)
        printf("I am proc 0 and the collected avg is %g, %g, %g\n", avg[0], avg[1], avg[2]);

    return 0;
}