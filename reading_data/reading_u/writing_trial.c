#include <stdio.h>
#include <string.h>
#include <netcdf.h>

// for write
#define FILE_NAME2 "map_summarized.nc"
#define UNITS_speed "m_s"
#define UNITS "units"
#define UNITS_time "s"
#define GRID_POINTS 8852366
#define TIME "time"  // the variable is called time in the file. I knew it using ncdump command
#define UNOD "unod"  // the variable is called unod in the file. I knew it using ncdump command
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); return 2;} //just simple macro to handle error with NETCDF
/*
mpicc -std=c99 -g -Wall -fopenmp -I /apps/netCDF4.7.0--gcc-9.1.0/include -L /apps/netCDF4.7.0--gcc-9.1.0/lib -lnetcdf -o writing_trial.out writing_trial.c -lm 

*/
int main (int argc, char *argv[]){
        /*Netcdf id for write*/
        int ncid2;
        int time_new_id;
        int gp_new_id;
        int time_var_new_id;
        int var_new_id;
        int y = 1;// indicating time step 1
        static float final_averages[GRID_POINTS] = {-1};
        int retval;
        if ((retval = nc_create(FILE_NAME2, NC_CLOBBER, &ncid2))) // ncclober to overwrite the file
        ERR(retval);
        /*define dimension*/
        if ((retval = nc_def_dim(ncid2, TIME, 1, &time_new_id)))
            ERR(retval);
        if ((retval = nc_def_dim(ncid2, UNOD, GRID_POINTS, &gp_new_id)))
            ERR(retval);
        /*defin variables*/
        if ((retval = nc_def_var(ncid2, "time", NC_INT, 1, &time_new_id, &time_var_new_id)))// define the variable
            ERR(retval);
        if ((retval = nc_def_var(ncid2, "speed", NC_FLOAT, 1, &gp_new_id, &var_new_id)))// define the variable
            ERR(retval);
        /*define text*/
        if ((retval = nc_put_att_text(ncid2, time_var_new_id, UNITS, 
				 strlen(UNITS_time), UNITS_time)))
            ERR(retval);
        if ((retval = nc_put_att_text(ncid2, var_new_id, UNITS, 
				 strlen(UNITS_speed), UNITS_speed)))
            ERR(retval);
        if ((retval = nc_enddef(ncid2)))
            ERR(retval);
        /*variables read*/
        if ((retval = nc_put_var_int(ncid2, time_var_new_id, &y)))
            ERR(retval);
        if ((retval = nc_put_var_float(ncid2, var_new_id, &final_averages[0])))
            ERR(retval);
        if ((retval = nc_close(ncid2)))
            ERR(retval);

}