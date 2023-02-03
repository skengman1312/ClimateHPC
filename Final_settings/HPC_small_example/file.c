/* the idea of this file is to create a small netcdf file to understand where is the problem is siutated */
/*
mpicc -g -Wall -I /apps/netCDF4.7.0--gcc-9.1.0/include -L /apps/netCDF4.7.0--gcc-9.1.0/lib -lnetcdf -o file.out file.c -lm -ldl -lz -lcurl -std=gnu99 -fopenmp
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <netcdf.h>
#include <mpi.h>
#include <omp.h>
#include <sys/time.h>
#define FILE_NAME2 "writing.nc"
#define TIME "time"  // the variable is called time in the file. I knew it using ncdump command
#define UNOD "unod"  // the variable is called unod in the file. I knew it using ncdump command
#define NZ "height"  // the variable is called unod in the file. I knew it using ncdump command

#define UNITS_speed "m_s"
#define UNITS "units"
#define UNITS_time "s"
#define GRID_POINTS 5
#define NZ1 10
#define ERR_spec(e) {printf("Error: %s\n", nc_strerror(e));} //just simple macro to handle error with NETCDF
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); return 2;} //just simple macro to handle error with NETCDF

#define NTIME 9

int main (int argc, char *argv[]){
    /*Creating 1 file for writing everything*/
    int retval,ncid2,time_dimid,speed_dimid,var_speed_id,level_dimid;
    int dimid[3]; // for wrtieing
    int i,j;
    /*START creating file */
    if ((retval = nc_create(FILE_NAME2, NC_CLOBBER, &ncid2))) // ncclober to overwrite the file
        ERR(retval);
    if ((retval = nc_def_dim(ncid2, TIME, NC_UNLIMITED, &time_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid2, NZ, NZ1, &level_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid2, UNOD, GRID_POINTS, &speed_dimid)))
        ERR(retval);
    dimid[0] = time_dimid;
    dimid[1] = level_dimid;
    dimid[2] = speed_dimid;
    if ((retval = nc_def_var(ncid2, "speed", NC_FLOAT, 3, dimid, &var_speed_id))) // define the varibael
        ERR(retval);
    if ((retval = nc_put_att_text(ncid2, var_speed_id, UNITS,
                                  strlen(UNITS_speed), UNITS_speed)))
        ERR(retval);
    /* End define mode. */
    if ((retval = nc_enddef(ncid2)))
        ERR(retval);
    /*Understanding how thing will proceed.*/
    float u_speed[NZ1][GRID_POINTS];
    for (i = 0; i < NZ1; i++)
    {        for (j = 0; j < GRID_POINTS; j++){
            u_speed[i][j] = i*GRID_POINTS+j;
            printf("%lf,",u_speed[i][j] );
        }
        printf("\n");
    }
    size_t start_1[3]={0,0,0};
    size_t count_1[3]={1,NZ1,GRID_POINTS};
    if ((retval = nc_put_vara_float(ncid2,var_speed_id, start_1, count_1,&u_speed[0][0])))
            ERR_spec(retval);

    for (i = 0; i < NZ1; i++)
    {        for (j = 0; j < GRID_POINTS; j++){
            u_speed[i][j] = i*GRID_POINTS+j+2;
            printf("%lf,",u_speed[i][j] );
        }
        printf("\n");
    }
    start_1[0] = 1;
    if ((retval = nc_put_vara_float(ncid2,var_speed_id, start_1, count_1,&u_speed[0][0])))
            ERR_spec(retval);

    for (i = 0; i < NZ1; i++)
    {        for (j = 0; j < GRID_POINTS; j++){
            u_speed[i][j] = i*GRID_POINTS+j+4;
            printf("%lf,",u_speed[i][j] );
        }
        printf("\n");
    }
    start_1[0] = 2;
    if ((retval = nc_put_vara_float(ncid2,var_speed_id, start_1, count_1,&u_speed[0][0])))
            ERR_spec(retval);


    for (i = 0; i < NZ1; i++)
    {        for (j = 0; j < GRID_POINTS; j++){
            u_speed[i][j] = i*GRID_POINTS+j+6;
            printf("%lf,",u_speed[i][j] );
        }
        printf("\n");
    }
    start_1[0] = 3;
    if ((retval = nc_put_vara_float(ncid2,var_speed_id, start_1, count_1,&u_speed[0][0])))
            ERR_spec(retval);

    if ((retval = nc_close(ncid2)))
        ERR(retval);
    /*END creating file */
}



