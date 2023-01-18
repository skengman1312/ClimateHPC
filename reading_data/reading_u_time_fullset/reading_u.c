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
/*FUNCTIONS START*/
/*This is a function that measures time using system time val
We subtract both time instances and we convert final output in seconds*/
double time_diff(struct timeval *start, struct timeval *end);
void convert_time_hour_sec( double seconds, long int *h, long int *m, long int *s);

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
    
    /*NETCDF id*/
    int ncid;
    int height_1_varid;          
    int time_varid;
    int unod_id;
    int retval;
    /*Netcdf id for write*/
    int ncid2;
    // int time_new_id;
    int gp_new_id;
    // int time_var_new_id;
    int var_new_id;
    int rec_dimid;
    double temp;
    /*Time variables to be used to see how much time each process takes*/
    struct timeval t_timer1_start;/*timer for process 0*/
    struct timeval t_timer1_finish;
    struct timeval t_timer2_start;
    struct timeval t_timer2_finish;
    struct timeval t_timer3_start;
    struct timeval t_timer3_finish;
    double t_nc_reading_time ;
    double t_threading_reading_time ;
    double t_nc_reading_time_sum ;
    double t_threading_reading_time_sum ;
    double t_nc_reading_time_Totalsum ;
    double t_threading_reading_time_Totalsum ;
    double t_comm_time ;
    double t_time_from_start;
    long int t_seconds = 0;
    long int t_minutes = 0;
    long int t_hours = 0;

  
    /* These program variables hold the time and depth. */
    double times[N_TIME], nz1s[N_NZ1]; 

    /* Program variables to hold the data we will read. 
    We will only need enough space to hold one timestep of data; one record. */
    // static float u_speed[GRID_POINTS]={-1};
    float *u_speed;
    u_speed = (float *)calloc(GRID_POINTS, sizeof(float));
    /*array of matrices containing all the different output average matrices calculated for each file*/
    // float *final_averages = malloc(N_TIME*GRID_POINTS* sizeof(float));
    float **final_averages = malloc(sizeof(float *) * N_TIME);
    for (int i = 0; i < N_TIME; i++)
        final_averages[i] = (float *)calloc(GRID_POINTS ,sizeof(float));
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
/*TIME 3 START*/
    if ((retval = nc_inq_varid(ncid, UNOD, &unod_id)))
            ERR(retval);
    if (rank == 0){
        gettimeofday(&t_timer3_start, NULL); //start timer of rank0
    }
    /*LOOOPING variables*/
    int rec;
    int i;
    int k ;
    /*Define the number of levels to do per process*/
    int levels_per_proc = ceil((double)N_NZ1 / size);/*if 2.3 is passed to ceil(), it will return 3*/
    /*you have 69level/4 proccess using ceil operator give you 18 level per process*/
    int limit = (rank + 1) * levels_per_proc;
    
    if (limit > N_NZ1)
    {
        limit = N_NZ1;
    }
    int thread_count;
#pragma omp parallel
    thread_count = omp_get_num_threads();
    /* Read the data. Since we know the contents of the file we know that the 
    data arrays in this program are the correct size to hold one timestep.*/
    count[0] = 1;/*1 time*/
    count[1] = 1;/*1 level*/
    count[2] = GRID_POINTS;/*all gridpoints*/
    // start[0] = k;
    start[1] = 0;
    start[2] = 0;

    /* end of setup of NetCDF reading */
    for (k = 0;k<N_TIME; k++){
        start[0] = k;
        t_nc_reading_time = 0 ;
        t_threading_reading_time = 0 ;
        t_nc_reading_time_sum = 0 ;
        t_threading_reading_time_sum = 0 ;
        t_nc_reading_time_Totalsum = 0;
        t_threading_reading_time_Totalsum = 0;
        t_comm_time = 0 ;
        t_time_from_start = 0;
        /* sum matrix */
        float * sum_u_speed;
        sum_u_speed = (float*)calloc(GRID_POINTS, sizeof(float));
        if (rank == 0)
        {
            printf("Number of processes: %d (levels being read for each process: %d)\n", size, levels_per_proc);
            printf("Number of threads: %d \n", thread_count);
/*TIME START T1*/
            gettimeofday(&t_timer1_start, NULL); //start timer of rank0
        }
        /* sum matrix */
        for (rec = rank * levels_per_proc; rec < limit; rec++){
/*TIME START T2*/
            gettimeofday(&t_timer2_start, NULL); //start reading timer
            start[1] = rec;
            if ((retval = nc_get_vara_float(ncid, unod_id, start, 
                        count, &u_speed[0])))
                ERR(retval);
            gettimeofday(&t_timer2_finish, NULL);
/*TIME END T2*/
            /* calculate the time take for reading*/
            t_nc_reading_time = time_diff(&t_timer2_start, &t_timer2_finish);
            t_nc_reading_time_sum+=t_nc_reading_time;
/*TIME START T2*/
            gettimeofday(&t_timer2_start, NULL); //start reading timer
    // #  pragma omp parallel for num_threads(thread_count)private(i )
            for (i = 0; i < GRID_POINTS;i++){
                sum_u_speed[i] += u_speed[i]/ N_NZ1;
            }
            gettimeofday(&t_timer2_finish, NULL);
/*TIME END T2*/
            t_threading_reading_time = time_diff(&t_timer2_start, &t_timer2_finish);
            t_threading_reading_time_sum+=t_threading_reading_time;
        }
    #ifdef DEBUG
        printf("The processes %d took %lf seconds to read all the nc data \n",rank,t_nc_reading_time_sum);
        printf("The processes %d took %lf seconds to thread \n",rank,t_threading_reading_time_sum);
    #endif
/*TIME START T2*/
    gettimeofday(&t_timer2_start, NULL); // start communication timer
    MPI_Reduce(&t_nc_reading_time_sum, &t_nc_reading_time_Totalsum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&t_threading_reading_time_sum, &t_threading_reading_time_Totalsum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(sum_u_speed, final_averages[k],GRID_POINTS, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);    
    gettimeofday(&t_timer2_finish, NULL);
/*TIME END T2*/
    t_comm_time=time_diff(&t_timer2_start, &t_timer2_finish);
    free(sum_u_speed);    
    if (rank == 0){ 
        gettimeofday(&t_timer1_finish, NULL); //start timer of rank0
/*TIME END T1*/
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
        printf("The end of readin time instance %d \n",k);
    }
    }

    if (rank == 0){
        gettimeofday(&t_timer3_finish, NULL);
/*TIME END 3*/

        temp=time_diff(&t_timer3_start, &t_timer3_finish);
        convert_time_hour_sec(temp,&t_hours,&t_minutes,&t_seconds);
        printf("The time taken to paralize everything for all of the files %lf seconds\n",temp);
        printf("The time taken to paralize everything for all of the files %ld hours,%ld minutes,%ld seconds \n",t_hours,t_minutes,t_seconds);
    }
    
    free(u_speed);    
    
    if (rank == 0)
    {
        /* Create the file. */
        // int y[12] = {1,2,3,4,5,6,7,8,9,10,11,12};                                            // indicating time step 1 
/*TIME START 2*/
        gettimeofday(&t_timer2_start, NULL); // start communication timer
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
        if ((retval = nc_put_vara_float(ncid2,var_new_id, start_1, count_1,&final_averages[rec][0])))
	            ERR(retval);
        }
        gettimeofday(&t_timer2_finish, NULL);
/*TIME END 2*/
        for (int i = 0; i < N_TIME; i++){
            free(final_averages[i]);      
        }
        if ((retval = nc_close(ncid2)))
            ERR(retval);
        temp=time_diff(&t_timer2_start, &t_timer2_finish);
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