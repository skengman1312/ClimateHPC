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
#define N_TIME 1            // the variable time has 12 entries
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
/*This is a function that measures time in minutes hours and seconds*/
void convert_time_hour_sec( double seconds, long int *h, long int *m, long int *s);

/*FUNCTIONS end*/

int main (int argc, char *argv[]){

    /* MPI  inizialization */
    MPI_Init(&argc, &argv);
    /*The total number of prcess that you will have*/
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    /*The rank of the current process*/
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    /*To save the processor name that I am working with*/
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    /*Print off a hello world message with processor name and rank*/
    printf("greetings:  %s, rank %d out of %d processes\n",
           processor_name, rank, size);

    /* VARIABLES DEFINE START */

   /*NETCDF id*/
    int ncid;
    int unod_id;
    int retval;
    
    /*Netcdf id for write*/
    int ncid2;
    int gp_new_id;
    int var_new_id;
    int rec_dimid;

    /*LOOOPING variables*/
    int rec;
    int i;
    int j;

    /* Variables related to OPmp*/
    int thread_count;
    int Time_per_proc;
    int limit;
    /*Time variables to be used to see how much time each process takes*/

    /*Stop watch 1 */
    struct timeval t_timer1_start;
    struct timeval t_timer1_finish;

    /*Stop watch 2 */
    struct timeval t_timer2_start;
    struct timeval t_timer2_finish;

    /*Used to output nc time taken for reach*/
    double t_nc_reading_time ;

    /*Used to output threading time of open mp*/
    double t_threading_reading_time ;

    /*Used to SUM output nc time taken for reach in 1 process*/
    double t_nc_reading_time_sum ;

    /*Used to SUM  output threading time of open mp in 1 process*/
    double t_threading_reading_time_sum ;

    /*Used to SUM output nc time taken for reach in ALL process*/
    double t_nc_reading_time_Totalsum ;

    /*Used to SUM  output threading time of open mp in ALL processes*/
    double t_threading_reading_time_Totalsum ;

    /*Time taken to do reduce operators over all process*/
    double t_comm_time ;

    /*time taken from beining of for loop until the end of the process.*/
    double t_time_from_start;

    /*To save number of seconds*/
    long int t_seconds = 0;

    /*To save number of minutes*/
    long int t_minutes = 0;

    /*To save number of hours*/
    long int t_hours = 0;
        
    /*For Netcdf reading proedure*/
    size_t start[NDIMS], count[NDIMS];

    /* VARIABLES DEFINE END*/

    /*DYNAMIC VARIABLES_intializaiton*/

    /*pointer of type float*/
    /*I want try to read full time stamp at once*/
    float **u_speed = malloc(sizeof(float *) * N_NZ1);
    // static float u_speed[N_NZ1][GRID_POINTS];
    // u_speed = (float *)calloc(GRID_POINTS, sizeof(float));
    for (int i = 0; i < N_NZ1; i++)
        u_speed[i] = (float *)calloc(GRID_POINTS ,sizeof(float));
    // static float final_averages[N_TIME][GRID_POINTS];

    float **final_averages = malloc(sizeof(float *) * N_TIME);
    for (int i = 0; i < N_TIME; i++)
        final_averages[i] = (float *)calloc(GRID_POINTS ,sizeof(float));
    
    float * sum_u_speed;
    sum_u_speed = (float*)calloc(GRID_POINTS, sizeof(float));
    /*DYNAMIC VARIABLES_intializaiton*/

    /* Open the file. */
    if ((retval = nc_open(FILE_NAME, NC_NOWRITE, &ncid)))
        ERR(retval);
    
    /* Get the varid of the velocity netCDF variable. */
    if ((retval = nc_inq_varid(ncid, UNOD, &unod_id)))
            ERR(retval);
    /*Define the number of levels to do per process*/
    Time_per_proc = ceil((double)N_TIME / size);/*if 2.3 is passed to ceil(), it will return 3*/
    /*you have 69level/4 proccess using ceil operator give you 18 level per process*/
    limit = (rank + 1) * Time_per_proc;
    if (limit > N_NZ1)
    {
        limit = N_NZ1;
    }
    /* Instialzing between the time variables */
    /*Get the number of threads*/
    #pragma omp parallel 
    thread_count = omp_get_num_threads();
    count[0] = 1;/*1 time*/
    count[1] = N_NZ1;/*ALL levels*/
    count[2] = GRID_POINTS;/*all gridpoints*/
    // start[0] = k;
    start[1] = 0;
    start[2] = 0;
    if (rank == 0)
        {
            printf("Number of processes: %d (Time being read for each process: %d)\n", size, Time_per_proc);
            printf("Number of threads: %d \n", thread_count);
            gettimeofday(&t_timer1_start, NULL); //start timer of rank0
        }
    for (rec = rank * Time_per_proc; rec < limit; rec++){
            gettimeofday(&t_timer2_start, NULL); //start reading timer
            start[0] = rec;
            if ((retval = nc_get_vara_float(ncid, unod_id, start, 
                        count, &u_speed[0][0])))
                ERR(retval);
            gettimeofday(&t_timer2_finish, NULL);
            /* calculate the time take for reading*/
            t_nc_reading_time = time_diff(&t_timer2_start, &t_timer2_finish);
            t_nc_reading_time_sum+=t_nc_reading_time;
            gettimeofday(&t_timer2_start, NULL); //start reading timer
    #  pragma omp parallel for num_threads(thread_count)private(i,j )
            for (i = 0; i < N_NZ1;i++){
                for (j = 0; j < GRID_POINTS;j++){
                    sum_u_speed[j] += u_speed[i][j]/ N_NZ1;
                }
            }
            gettimeofday(&t_timer2_finish, NULL);
            t_threading_reading_time = time_diff(&t_timer2_start, &t_timer2_finish);
            t_threading_reading_time_sum+=t_threading_reading_time;
        }
    #ifdef DEBUG
        printf("The processes %d took %lf seconds to read all the nc data \n",rank,t_nc_reading_time_sum);
        printf("The processes %d took %lf seconds to thread \n",rank,t_threading_reading_time_sum);
    #endif
    gettimeofday(&t_timer2_start, NULL); // start communication timer
    MPI_Reduce(&t_nc_reading_time_sum, &t_nc_reading_time_Totalsum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&t_threading_reading_time_sum, &t_threading_reading_time_Totalsum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(sum_u_speed, final_averages,GRID_POINTS * N_TIME, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);    
    gettimeofday(&t_timer2_finish, NULL);
    t_comm_time=time_diff(&t_timer2_start, &t_timer2_finish);
    free(sum_u_speed);    

    if (rank == 0){ 
        gettimeofday(&t_timer1_finish, NULL); //start timer of rank0
        t_time_from_start=time_diff(&t_timer1_start, &t_timer1_finish);
        convert_time_hour_sec(t_nc_reading_time_Totalsum,&t_hours,&t_minutes,&t_seconds);
        printf("The time taken to do Nc read is %lf seconds\n",t_nc_reading_time_Totalsum);
        printf("The time taken to do Nc read is %ld hours,%ld minutes,%ld seconds \n", t_hours,t_minutes,t_seconds);
        /**/
        convert_time_hour_sec(t_threading_reading_time_Totalsum,&t_hours,&t_minutes,&t_seconds);
        printf("The time taken to do the threading is %lf seconds\n",t_threading_reading_time_Totalsum);
        printf("The time taken to do the threading is %ld hours,%ld minutes,%ld seconds \n",t_hours,t_minutes,t_seconds);
        /**/
        convert_time_hour_sec(t_comm_time,&t_hours,&t_minutes,&t_seconds);
        printf("The time taken to do 3 MPI REDUCE is %lfseconds \n",t_comm_time);
        printf("The time taken to do 3 MPI REDUCE is %ld hours,%ld minutes,%ld seconds \n",t_hours,t_minutes,t_seconds);
        /**/
        convert_time_hour_sec(t_time_from_start,&t_hours,&t_minutes,&t_seconds);
        printf("The time taken from start of For loop till the reduce is %lf seconds\n",t_time_from_start);
        printf("The time taken from start of For loop till the reduce is %ld hours,%ld minutes,%ld seconds \n",t_hours,t_minutes,t_seconds);
    }
    for (int i = 0; i < N_TIME; i++){
            free(u_speed[i]);      
    }
    if (rank == 0)
    {
        /* Create the file. */
        if ((retval = nc_create(FILE_NAME2, NC_CLOBBER, &ncid2))) // ncclober to overwrite the file
            ERR(retval);
        /*define dimension*/
        // if ((retval = nc_def_dim(ncid2, TIME, 12, &time_new_id)))
        //     ERR(retval);
        if ((retval = nc_def_dim(ncid2, TIME, NC_UNLIMITED, &rec_dimid)))
            ERR(retval);
        if ((retval = nc_def_dim(ncid2, UNOD, GRID_POINTS, &gp_new_id)))
            ERR(retval);
        int dimid[2];
        dimid[0]=rec_dimid;
        dimid[1]=gp_new_id;
        /*define variable*/
        if ((retval = nc_def_var(ncid2, "speed", NC_FLOAT, 2, dimid, &var_new_id)))// define the varibael
            ERR(retval);
        // if ((retval = nc_def_var(ncid2, "time", NC_INT, 1, &time_new_id, &time_var_new_id)))// define the varibael
        //     ERR(retval);
        // if ((retval = nc_put_att_text(ncid2, time_var_new_id, UNITS, 
		// 		 strlen(UNITS_time), UNITS_time)))
        //     ERR(retval);
        if ((retval = nc_put_att_text(ncid2, var_new_id, UNITS, 
				 strlen(UNITS_speed), UNITS_speed)))
            ERR(retval);
        /* End define mode. */
        if ((retval = nc_enddef(ncid2)))
            ERR(retval);
        // if ((retval = nc_put_var_int(ncid2, time_var_new_id, &y[0])))
        //     ERR(retval);
        size_t start_1[2];
        size_t count_1[2];
        count_1[0] = 1;
        count_1[1] = GRID_POINTS;
        start_1[1] = 0;

        for (rec = 0; rec < N_TIME; rec++)
        {
        start_1[0] = rec;
        if ((retval = nc_put_vara_float(ncid2,var_new_id, start_1, count_1,&final_averages[0][0])))
	            ERR(retval);
        }
        for (int i = 0; i < N_TIME; i++){
            free(final_averages[i]);      
        }
        if ((retval = nc_close(ncid2)))
            ERR(retval);
        gettimeofday(&t_timer2_finish, NULL);
        double temp=time_diff(&t_timer2_start, &t_timer2_finish);
        convert_time_hour_sec(temp,&t_hours,&t_minutes,&t_seconds);
        printf("The time taken to wrtie the file is %lf seconds\n",temp);
        printf("The time taken to wrtie the file is %ld hours,%ld minutes,%ld seconds \n",t_hours,t_minutes,t_seconds);
    }
    /*CLOSING FILE*/
    if ((retval = nc_close(ncid)))
        ERR(retval);
    printf("This process has about to close right:on this node %s,with this ranking %d out of %d processes \n",
           processor_name, rank, size);
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