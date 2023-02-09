#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <netcdf.h>
#include <mpi.h>
#include <omp.h>
#include <sys/time.h>


/*MACROS START*/
#define FILE_NAME "/shares/HPC4DataScience/FESOM2/unod.fesom.2010.nc"
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); return 2;} //just simple macro to handle error with NETCDF
#define ERR_spec(e) {printf("Error: %s\n", nc_strerror(e));} //just simple macro to handle error with NETCDF
#define HEIGHT_1_NAME "nz1"  // the variable is called nz in the file. I knew it using ncdump command
#define TIME "time"  // the variable is called time in the file. I knew it using ncdump command
#define UNOD "unod"  // the variable is called unod in the file. I knew it using ncdump command
#define N_NZ1 69            // the variable nz1 has 69 entries
#define N_TIME 12            // the variable time has 12 entries
#define NDIMS 3
#define GRID_POINTS 8852366
#define DEBUG 1
// for write
#define FILE_NAME2 "map_summarized_fuldataset_serial.nc"
#define UNITS_speed "m_s"
#define UNITS "units"
#define UNITS_time "s"
#define NDIMS_wr 3
/*MACROS end*/

double time_diff(struct timeval *start, struct timeval *end);
void convert_time_hour_sec( double seconds, long int *h, long int *m, long int *s);
void net_write(float * final_averages,int k);

int main (int argc, char *argv[]){

    /*LOOOPING variables*/
    int rec;
    int i;
    int k ;
    /*NETCDF id*/
    int ncid;
    int unod_id;
    int retval;
    /*Netcdf id for write*/
    int ncid2; // write file nmae
    int time_dimid; // time dimension id 
    int speed_dimid; // speed dimension id 
    int dimid[2]; // for wrtieing 
    int var_speed_id;// variable speed id 
    size_t start[NDIMS], count[NDIMS];

    /*Time variables to be used to see how much time each process takes*/
    struct timeval t_timer1_start;/*timer for process 0*/
    struct timeval t_timer1_finish;
    struct timeval t_timer3_start;
    struct timeval t_timer3_finish;
    double t_time_from_start;
    long int t_seconds = 0;
    long int t_minutes = 0;
    long int t_hours = 0;
    double temp;
    float *u_speed;
    u_speed = (float *)calloc(GRID_POINTS, sizeof(float));
    /*Creating 1 file for writing everything*/
        /*START creating file */
        if ((retval = nc_create(FILE_NAME2, NC_CLOBBER, &ncid2))) // ncclober to overwrite the file
            ERR(retval);
        if ((retval = nc_def_dim(ncid2, TIME, NC_UNLIMITED, &time_dimid)))
            ERR(retval);
        if ((retval = nc_def_dim(ncid2, UNOD, GRID_POINTS, &speed_dimid)))
            ERR(retval);
        dimid[0]=time_dimid;
        dimid[1]=speed_dimid;
        if ((retval = nc_def_var(ncid2, "speed", NC_FLOAT, 2, dimid, &var_speed_id)))// define the varibael
            ERR(retval);
        if ((retval = nc_put_att_text(ncid2, var_speed_id, UNITS, 
				 strlen(UNITS_speed), UNITS_speed)))
            ERR(retval);
        /* End define mode. */
        if ((retval = nc_enddef(ncid2)))
            ERR(retval);
        if ((retval = nc_close(ncid2)))
            ERR(retval);
        /*END creating file */
    
    if ((retval = nc_open(FILE_NAME, NC_NOWRITE, &ncid)))
    ERR(retval);
    if ((retval = nc_inq_varid(ncid, UNOD, &unod_id)))
    ERR(retval);
    count[0] = 1;/*1 time*/
    count[1] = 1;/*1 level*/
    count[2] = GRID_POINTS;/*all gridpoints*/
    start[1] = 0;
    start[2] = 0;
    /* START timer 3 to have calculate total time*/
    gettimeofday(&t_timer3_start, NULL); //start timer of rank0
    for (k = 0;k<N_TIME; k++){

        start[0] = k;
        t_time_from_start = 0;
        float * sum_u_speed;
        sum_u_speed = (float*)calloc(GRID_POINTS, sizeof(float));
        gettimeofday(&t_timer1_start, NULL); //start timer of rank0
        for (rec = 0; rec < N_NZ1; rec++){
            start[1] = rec;
            if ((retval = nc_get_vara_float(ncid, unod_id, start, count, &u_speed[0])))
                ERR(retval);
            for (i = 0; i < GRID_POINTS;i++){
                sum_u_speed[i] += u_speed[i]/ N_NZ1;
            }
        }
        gettimeofday(&t_timer1_finish, NULL); //start timer of rank0
        /*TIME END T1*/
        t_time_from_start=time_diff(&t_timer1_start, &t_timer1_finish);
        /**/
        convert_time_hour_sec(t_time_from_start,&t_hours,&t_minutes,&t_seconds);
        printf("##### THE BEGINING OF THE RESULT OF INSTANCE %d ##### \n",k);
        printf("The time taken from start of For loop till the reduce is %lf seconds\n", t_time_from_start);
        printf("The time taken from start of For loop till the reduce is %ld hours,%ld minutes,%ld seconds \n",t_hours,t_minutes,t_seconds);
        printf("##### THE END OF THE RESULT OF INSTANCE %d #####\n ",k);
        printf("\n");
        net_write(sum_u_speed, k);
        free(sum_u_speed);

    }
    gettimeofday(&t_timer3_finish, NULL);
    /*TIME END 3*/
    temp=time_diff(&t_timer3_start, &t_timer3_finish);
    convert_time_hour_sec(temp,&t_hours,&t_minutes,&t_seconds);
    printf("####THIS IS THE FINAL RESULTSSSS###\n");
    printf("The time taken to paralize everything for all of the files %lf seconds\n",temp);
    printf("The time taken to paralize everything for all of the files %ld hours,%ld minutes,%ld seconds \n",t_hours,t_minutes,t_seconds);
    /*CLOSING FILE*/
    if ((retval = nc_close(ncid)))
        ERR(retval);
    return 0;
}

void net_write(float * final_averages, int k){
   int ncid, retval,unod_id;
    size_t start_1[2]={k,0};
    size_t count_1[2]={1,GRID_POINTS};
    if ((retval = nc_open(FILE_NAME2, NC_WRITE, &ncid)))ERR_spec(retval);
    if ((retval = nc_inq_varid(ncid, "speed", &unod_id)))ERR_spec(retval);
    if ((retval = nc_put_vara_float(ncid,unod_id, start_1, count_1,&final_averages[0])))ERR_spec(retval);
    if ((retval = nc_close(ncid)))ERR_spec(retval);
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