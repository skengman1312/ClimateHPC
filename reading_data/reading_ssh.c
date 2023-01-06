#include <stdio.h>
#include <string.h>
#include <netcdf.h>
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
/*MACROS end*/
int main (){
    // /* VARIABLES DEFINE START*/
    /*NETCDF id*/
    int ncid;
    int ssh_varid;
    int time_varid;
    int ssh_id;
    int retval;

    /* These program variables hold the time and depth. */
    double times[N_TIME], nz1s[N_NZ1];

    /* Program variables to hold the data we will read. 
    We will need enough space to hold all the timesteps of data. */
    static float ssh[N_TIME][GRID_POINTS];

    /* The start and count arrays will tell the netCDF library where to read our data. */
    size_t start[NDIMS], count[NDIMS];

    /* Open the file. */
    if ((retval = nc_open(FILE_NAME, NC_NOWRITE, &ncid)))
        ERR(retval);

    /* Get the varids of the time and nz1 coordinate variables. */
    if ((retval = nc_inq_varid(ncid, SSH, &ssh_varid)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, TIME, &time_varid)))
        ERR(retval);


    /* Get the varid of the ssh netCDF variable. */
    if ((retval = nc_inq_varid(ncid, SSH, &ssh_id)))
            ERR(retval);
    /* Read the data. Since we know the contents of the file we know that the 
    data arrays in this program are the correct size to hold one timestep.*/
    count[0] = 1;
//    count[1] = N_NZ1;
    count[1] = GRID_POINTS;
    start[0] = 0;
    start[1] = 0;
    /* end of setup of NetCDF reading */
    
    /*LOOOPING variables*/
    int rec; 
    /* sum matrix */
//    static float sum_u_speed[N_NZ1][GRID_POINTS] = {{0}};
    for (rec = 0; rec < 30; rec++){
        start[0] = rec;
        if ((retval = nc_get_vara_float(ncid, ssh_id, start,
				      count, &ssh[rec][0])))
            ERR(retval);

    }
    printf("%lf\n",ssh[0][0]);
    printf("%lf\n",ssh[29][8852365]);
    /*CLOSING FILE*/
    if ((retval = nc_close(ncid)))
        ERR(retval);
    printf("*** SUCCESS :) reading example %s\n", FILE_NAME);


    return 0;
}

