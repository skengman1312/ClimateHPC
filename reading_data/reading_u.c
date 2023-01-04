#include <stdio.h>
#include <string.h>
#include <netcdf.h>
#include <mpi.h>
#include <omp.h>

/*
mpicc -std=c99 -g -Wall -I /apps/netCDF4.7.0--gcc-9.1.0/include -L /apps/netCDF4.7.0--gcc-9.1.0/lib -lnetcdf -o reading_u.out reading_u.c -lm 
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
/*MACROS end*/
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
    printf("greetings:  %s, rank %d out of %d processors\n",
           processor_name, rank, size);
    /* VARIABLES DEFINE START*/
    
    /*NETCDF id*/
    int ncid;
    int height_1_varid;          
    int time_varid;
    int unod_id;
    int retval;

    /* These program variables hold the time and depth. */
    double times[N_TIME], nz1s[N_NZ1]; 

    /* Program variables to hold the data we will read. 
    We will only need enough space to hold one timestep of data; one record. */
    static float u_speed[GRID_POINTS];

    /* The start and count arrays will tell the netCDF library where to read our data. */
    size_t start[NDIMS], count[NDIMS];
    /* VARIABLES DEFINE END*/

    /* Open the file. */
    if ((retval = nc_open(FILE_NAME, NC_NOWRITE, &ncid)))
        ERR(retval);

    /* Get the varids of the time and nz1 coordinate variables. */
    if ((retval = nc_inq_varid(ncid, HEIGHT_1_NAME, &height_1_varid)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, TIME, &time_varid)))
        ERR(retval);

    /* Read the coordinate variable data. */
    if ((retval = nc_get_var_double(ncid, height_1_varid, &nz1s[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, time_varid, &times[0])))
        ERR(retval);

    /* Get the varid of the velocity netCDF variable. */
    if ((retval = nc_inq_varid(ncid, UNOD, &unod_id)))
            ERR(retval);
    /* Read the data. Since we know the contents of the file we know that the 
    data arrays in this program are the correct size to hold one timestep.*/
    count[0] = 1;
    count[1] = 1;
    count[2] = GRID_POINTS;
    start[0] = 0;
    // start[1] = 0;
    start[2] = 0;
    /* end of setup of NetCDF reading */
    
    /*LOOOPING variables*/
    int rec; 
    /* sum matrix */
    // static float sum_u_speed[N_NZ1][GRID_POINTS] = {{0}};
    for (rec = 0; rec < N_NZ1; rec++){
        start[1] = rec;
        if ((retval = nc_get_vara_float(ncid, unod_id, start, 
				      count, &u_speed[0])))
            ERR(retval);
        printf("%lf\n",u_speed[0]);
        break;
    }
    /*CLOSING FILE*/
    if ((retval = nc_close(ncid)))
        ERR(retval);
    printf("*** SUCCESS reading example unod.fesom.2010.nc!\n");
    MPI_Finalize();
    return 0;
}

