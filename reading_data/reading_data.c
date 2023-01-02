#include <stdio.h>
#include <string.h>
#include <netcdf.h>
/*
mpicc -std=c99 -g -Wall -I /apps/netCDF4.7.0--gcc-9.1.0/include -L /apps/netCDF4.7.0--gcc-9.1.0/lib -lnetcdf -o reading_data_1 reading_data.c -lm 
*/
/*MACROS START*/
#define FILE_NAME "/shares/HPC4DataScience/FESOM2/fesom.mesh.diag.nc"
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); return 2;} //just simple macro to handle error with NETCDF
#define HEIGHT_1_NAME "nz"  // the variable is called nz in the file. I knew it using ncdump command
#define HEIGHT_2_NAME "nz1" // the variable is called nz1 in the file. I knew it using ncdump command
#define LONGITUDE_NAME "lon"
#define LATITUDE_NAME "lat"
#define N_NZ 70             // the variable nz has 70 entries
#define N_NZ1 69            // the variable nz1 has 69 entries
#define START_nz 0
#define START_nz1 -2.5
#define NDIMS 3 // Variable to tell netcdf how many dimenisons we have
#define GRID_POINTS 8852366
/*MACROS end*/
int main (){
    /* VARIABLES DEFINE START*/
    int ncid;               // nc file id  
    int height_1_varid;          // nz id
    int height_2_varid;          // nz1 id
    int lat_varid;              // Lat id
    int lon_varid;              // long id 
    int retval;             // return value to be printed
    double nzs[N_NZ], nz1s[N_NZ1]; // two arrays .
    /* The start and count arrays will tell the netCDF library where to
      read our data. */
    size_t start[NDIMS], count[NDIMS];
    /* VARIABLES DEFINE END*/
    /* Open the file. */
    if ((retval = nc_open(FILE_NAME, NC_NOWRITE, &ncid)))
        ERR(retval);
    /* Get the varids of the nz and nz1 coordinate
    * variables. */
    if ((retval = nc_inq_varid(ncid, HEIGHT_1_NAME, &height_1_varid)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, HEIGHT_2_NAME, &height_2_varid)))
        ERR(retval);
    // // /* Read the coordinate variable data. */
    if ((retval = nc_get_var_double(ncid, height_1_varid, &nzs[0])))
        ERR(retval);
    if ((retval = nc_get_var_double(ncid, height_2_varid, &nz1s[0])))
        ERR(retval);
    /* Get the varids of the LONG and LATITUDE netCDF
    * variables. */
    if ((retval = nc_inq_varid(ncid, LONGITUDE_NAME, &lon_varid)))
        ERR(retval);
    if ((retval = nc_inq_varid(ncid, LATITUDE_NAME, &lat_varid)))
        ERR(retval);
    count[0] = N_NZ;
    count[1] = N_NZ1;
    count[2] = GRID_POINTS;
    start[0] = 0;
    start[1] = 0;
    start[2] = 0;
    // if ((retval = nc_get_vara_double(ncid, pres_varid, start, count, &pres_in[0][0][0])))
	//     ERR(retval);
    /*CLOSING FILE*/
    if ((retval = nc_close(ncid)))
        ERR(retval);
    printf("*** SUCCESS reading example file pres_temp_4D.nc!\n");
    return 0;
}

