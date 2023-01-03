#include <stdio.h>
#include <string.h>
#include <netcdf.h>
/*
mpicc -std=c99 -g -Wall -I /apps/netCDF4.7.0--gcc-9.1.0/include -L /apps/netCDF4.7.0--gcc-9.1.0/lib -lnetcdf -o reading_data reading_data.c -lm 
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
int main (){
    // /* VARIABLES DEFINE START*/
    int ncid;               // nc file id  
    int height_1_varid;          // nz1 id
    int time_varid;             // time id
    int unod_id;                        //unod id
    int retval; // return value to be printed
    int rec; // variable for looping
    double times[N_TIME], nz1s[N_NZ1]; // two arrays containing the index .
    static float u_speed[N_NZ1][GRID_POINTS];
    // /* The start and count arrays will tell the netCDF library where to
    //   read our data. */
    size_t start[NDIMS], count[NDIMS];
    // /* VARIABLES DEFINE END*/
    // /* Open the file. */
    if ((retval = nc_open(FILE_NAME, NC_NOWRITE, &ncid)))
        ERR(retval);
    // /* Get the varids of the nz and nz1 coordinate
    // * variables. */
    if ((retval = nc_inq_varid(ncid, HEIGHT_1_NAME, &height_1_varid)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, TIME, &time_varid)))
        ERR(retval);
    // // // /* Read the coordinate variable data. */
    if ((retval = nc_get_var_double(ncid, height_1_varid, &nz1s[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, time_varid, &times[0])))
        ERR(retval);
    // This function is only used for testing delet it later START
    // int i;
    // for (i = 0; i < 12; i++)
    // {
    //     printf("%lf,",times[i] );
    // }
    // This function is only used for testing delet it later END

    // /* Get the varids of the velocity
    // * variables. */
    if ((retval = nc_inq_varid(ncid, UNOD, &unod_id)))
            ERR(retval);
    count[0] = 1;
    count[1] = N_NZ1;
    count[2] = GRID_POINTS;
    start[1] = 0;
    start[2] = 0;
    for (rec = 0; rec < N_TIME; rec++){
        start[0] = rec;
        if ((retval = nc_get_vara_float(ncid, unod_id, start, 
				      count, &u_speed[0][0])))
            ERR(retval);
        printf("%lf",u_speed[0][0]);
        break;
    }
    /*CLOSING FILE*/
    if ((retval = nc_close(ncid)))
        ERR(retval);
    printf("*** SUCCESS reading example unod.fesom.2010.nc!\n");
    return 0;
}

