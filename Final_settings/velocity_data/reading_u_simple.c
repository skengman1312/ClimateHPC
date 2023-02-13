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
#define FILE_NAME2 "map_summarized_full_dataset_simple_version.nc"
#define UNITS_speed "m_s"
#define UNITS "units"
#define UNITS_time "s"
#define NDIMS_wr 3
/*MACROS end*/

double time_diff(struct timeval *start, struct timeval *end);
void convert_time_hour_sec( double seconds, long int *h, long int *m, long int *s);
void net_write(float * final_averages,int k);

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
    printf("greetings:  %s, rank %d out of %d processes\n",processor_name, rank, size);
    /* VARIABLES DEFINE START*/
    
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
    struct timeval t_timer2_start;
    struct timeval t_timer2_finish;
    struct timeval t_timer3_start;
    struct timeval t_timer3_finish;
    struct timeval starttime;
    struct timeval endtime;

    double walltimes_start[2];
    double walltimes_end[2];
    double t_comm_time ;
    double t_time_from_start;
    long int t_seconds = 0;
    long int t_minutes = 0;
    long int t_hours = 0;
    double temp;
    double passing_time = 0;
    double nc_reading ; //variable containing the local sum of all elapsed time for reading
    double threading_time = 0;  //variable containing the local sum of all elapsed time for writing the sum matrix
    double sum_nc_reading;        // variable containing the output of reduce operation, collecting all reading times of different processes
    double sum_threading_time; // variable containing the output of reduce operation, collecting all elaboration times of different processes
    double total_time_reading=0;
    double total_time_threading=0;
    float *u_speed;
    u_speed = (float *)calloc(GRID_POINTS, sizeof(float));
    /*Creating 1 file for writing everything*/
    if(0==rank){
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
    }

    if ((retval = nc_open(FILE_NAME, NC_NOWRITE, &ncid)))
    ERR(retval);
    if ((retval = nc_inq_varid(ncid, UNOD, &unod_id)))
    ERR(retval);
    int levels_per_proc = ceil((double)N_NZ1 / size);/*if 2.3 is passed to ceil(), it will return 3*/

    int limit = (rank + 1) * levels_per_proc;
    if (limit > N_NZ1)
    {
        limit = N_NZ1;
    }
    count[0] = 1;/*1 time*/
    count[1] = 1;/*1 level*/
    count[2] = GRID_POINTS;/*all gridpoints*/
    start[1] = 0;
    start[2] = 0;
    /* START timer 3 to have calculate total time*/
    if (rank == 0){
        gettimeofday(&t_timer3_start, NULL); //start timer of rank0
        walltimes_start[1] = MPI_Wtime();
    }

    for (k = 0;k<N_TIME; k++){
        start[0] = k;
        t_comm_time = 0 ;
        t_time_from_start = 0;
        nc_reading = 0;
        threading_time = 0;
        float *sum_u_speed;
        sum_u_speed = (float*)calloc(GRID_POINTS, sizeof(float));
        float *final_averages;
        final_averages = (float *)calloc(GRID_POINTS, sizeof(float));
        if (rank == 0){
            printf("Number of processes: %d (levels being read for each process: %d)\n", size, levels_per_proc);
            /*TIME START T1*/
            gettimeofday(&t_timer1_start, NULL); //start timer of rank0
            walltimes_start[0] = MPI_Wtime();

        }
        for (rec = rank * levels_per_proc; rec < limit; rec++){
            start[1] = rec;
            gettimeofday(&starttime, NULL); //start reading timer
            if ((retval = nc_get_vara_float(ncid, unod_id, start, count, &u_speed[0])))
                ERR(retval);
            gettimeofday(&endtime, NULL); //start reading timer
            passing_time = time_diff(&starttime, &endtime);
            nc_reading += passing_time;
            gettimeofday(&starttime, NULL); //start reading timer
            for (i = 0; i < GRID_POINTS;i++){
                sum_u_speed[i] += u_speed[i]/ N_NZ1;
            }
            gettimeofday(&endtime, NULL); //start reading timer
            passing_time = time_diff(&starttime, &endtime);
            threading_time += passing_time;
        }
        MPI_Reduce(&nc_reading, &sum_nc_reading, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&threading_time, &sum_threading_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        /*TIME START T2*/
        gettimeofday(&t_timer2_start, NULL); // start communication timer
        MPI_Reduce(sum_u_speed, final_averages,GRID_POINTS, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);    
        gettimeofday(&t_timer2_finish, NULL);
        /*TIME END T2*/

        /*Communication time*/
        t_comm_time=time_diff(&t_timer2_start, &t_timer2_finish);

        if (rank == 0){ 
            gettimeofday(&t_timer1_finish, NULL); //start timer of rank0
            walltimes_end[0] = MPI_Wtime();
            /*TIME END T1*/
            net_write(final_averages,k);
            total_time_reading += sum_nc_reading / size;
            total_time_threading += sum_threading_time / size;
            t_time_from_start=time_diff(&t_timer1_start, &t_timer1_finish);
            /*Total time i did the communication between the process*/
            printf("##### THE BEGINING OF THE RESULT OF INSTANCE %d ##### \n",k);
            convert_time_hour_sec(t_comm_time,&t_hours,&t_minutes,&t_seconds);
            printf("The time taken to do 3 MPI REDUCE is %lfseconds \n",t_comm_time);
            printf("The time taken to do 3 MPI REDUCE is %ld hours,%ld minutes,%ld seconds \n",t_hours,t_minutes,t_seconds);
            
            printf("Average reading time per process: %.7f\n", sum_nc_reading / size);
            printf("Average threading time per process: %.7f\n\n", sum_threading_time / size);
            /*Total time for the code*/
            convert_time_hour_sec(t_time_from_start,&t_hours,&t_minutes,&t_seconds);
            printf("The time taken from start of For loop till the reduce for 1 time step is %lf seconds\n",t_time_from_start);
            printf("The time taken from start of For loop till the reduce for 1 time step is %ld hours,%ld minutes,%ld seconds \n",t_hours,t_minutes,t_seconds);
            printf("Total walltime from the begining till the end for 1 time step is %lf seconds \n", walltimes_end[0] - walltimes_start[0]);
            printf("##### THE END OF THE RESULT OF INSTANCE %d #####\n ",k);
            printf("\n");        
        }
        free(sum_u_speed);
        free(final_averages);
    }
    if (rank == 0){
        gettimeofday(&t_timer3_finish, NULL);
        walltimes_end[1] = MPI_Wtime();
        /*TIME END 3*/
        temp=time_diff(&t_timer3_start, &t_timer3_finish);
        convert_time_hour_sec(temp,&t_hours,&t_minutes,&t_seconds);
        printf("####THIS IS THE FINAL RESULTSSSS###\n");
        printf("Average reading time per process over 12 month: %.7f\n", total_time_reading / N_TIME);
        printf("Average threading time per process over 12 month: %.7f\n\n", total_time_threading / N_TIME);
        printf("The time taken to paralize everything for all of the files %lf seconds\n",temp);
        printf("The time taken to paralize everything for all of the files %ld hours,%ld minutes,%ld seconds \n",t_hours,t_minutes,t_seconds);
        printf("Total walltime from the begining till the end for ALL time step is %lf seconds \n", walltimes_end[1] - walltimes_start[1]);
    }
    free(u_speed);
    /*CLOSING FILE*/
    if ((retval = nc_close(ncid)))
        ERR(retval);
    printf("This process has about to close right:on this node %s,with this ranking %d out of %d processes \n",
           processor_name, rank, size);
    MPI_Finalize();
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