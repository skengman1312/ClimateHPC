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
#define FILE_NAME2 "map_summarized_fuldataset.nc"
#define UNITS_speed "m_s"
#define UNITS "units"
#define UNITS_time "s"
#define NDIMS_wr 3
/*MACROS end*/
/*FUNCTIONS START*/
/*This is a function that measures time using system time val
We subtract both time instances and we convert final output in seconds*/
double time_diff(struct timeval * start, struct timeval *end);
void convert_time_hour_sec( double seconds, long int *h, long int *m, long int *s);
void threading(float * sum, float levels[][GRID_POINTS], int rank, int size, int levels_per_proc, int thread_count);
void net_write(float * final_averages,int k);
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
    /*NETCDF id*/
    int ncid;
    int unod_id;
    int retval;
    /*Netcdf id for write*/
    int ncid2; // write file nmae
    int time_dimid; // time dimension id 
    int speed_dimid; // speed dimension id 
    int dimid[2];
    int var_speed_id;// variable speed id 
    size_t start[NDIMS], count[NDIMS];
    /*LOOOPING variables*/
    // int rec;
    int i;
    int k;
    int levels_per_proc = ceil((double)N_NZ1 / size);/*if 2.3 is passed to ceil(), it will return 3*/
    int sendcounts[size];// for scatterv
    int displs[size];// for scatterv
    int thread_count;
    double temp;
    double temp2;
    // intializing the arrays
    for (i = 0; i < size;i++){
        sendcounts[i]= levels_per_proc * GRID_POINTS;
        displs[i] = i * levels_per_proc * GRID_POINTS;
    }
    sendcounts[size-1] = (N_NZ1-(levels_per_proc*(size-1)))*GRID_POINTS;
    static float u_speed[N_NZ1][GRID_POINTS] = {{0}};
    float levels[levels_per_proc][GRID_POINTS];
    #pragma omp parallel
    thread_count = omp_get_num_threads();
    /*Time variables to be used to see how much time each process takes*/
    struct timeval t_timer1_start;/*timer for process 0*/
    struct timeval t_timer1_finish;
    struct timeval t_timer2_start;
    struct timeval t_timer2_finish;
    struct timeval t_timer3_start;
    struct timeval t_timer3_finish;

    long int t_seconds = 0;
    long int t_minutes = 0;
    long int t_hours = 0;
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
        if ((retval = nc_open(FILE_NAME, NC_NOWRITE, &ncid)))
        ERR(retval);
        if ((retval = nc_inq_varid(ncid, UNOD, &unod_id)))
        ERR(retval);
        count[0] = 1;/*1 time*/
        count[1] = N_NZ1;/*all level*/
        count[2] = GRID_POINTS;/*all gridpoints*/
        start[1] = 0;
        start[2] = 0;
    }
    /*START*/
    if (rank == 0)
        {
            printf("Number of processes: %d (levels being read for each process: %d)\n", size, levels_per_proc);
            printf("Number of threads: %d \n", thread_count);
/*TIME START T3*/
            printf("#########THE START OF COMPUTATION OVER ALL TIMESTEPS#######\n");
            gettimeofday(&t_timer3_start, NULL); //start timer of rank0
        }
    for (k = 0;k<N_TIME; k++){
        float * sum;
        float *final_averages;
        sum = (float *)calloc(GRID_POINTS, sizeof(float));
        final_averages = (float *)calloc(GRID_POINTS, sizeof(float));
        if (rank == 0)
        {
            gettimeofday(&t_timer1_start, NULL); //start timer of rank0

        }
        if(rank==0){
            start[0] = k;
            if ((retval = nc_get_vara_float(ncid, unod_id, start, count, &u_speed[0][0])))
                        ERR(retval);
            /*SCATTERING*/
            gettimeofday(&t_timer2_start, NULL); // start communication timer
            MPI_Scatterv(u_speed,sendcounts,displs,MPI_FLOAT,levels,GRID_POINTS*levels_per_proc,MPI_FLOAT,0,MPI_COMM_WORLD);
            gettimeofday(&t_timer2_finish, NULL);
            temp=time_diff(&t_timer2_start, &t_timer2_finish);
            #ifdef DEBUG
            convert_time_hour_sec(temp,&t_hours,&t_minutes,&t_seconds);
            printf("The processes %d took %lf seconds to scatter \n",rank,temp);
            printf("The process took this time to finish scattering %ld hours,%ld minutes,%ld seconds \n",t_hours,t_minutes,t_seconds);
            #endif
            /*THREADING*/
            gettimeofday(&t_timer2_start, NULL); // start communication timer
            threading(sum, levels, rank, size, levels_per_proc, thread_count);
            gettimeofday(&t_timer2_finish, NULL);
            temp=time_diff(&t_timer2_start, &t_timer2_finish);
            #ifdef DEBUG
            convert_time_hour_sec(temp,&t_hours,&t_minutes,&t_seconds);
            printf("The processes %d took %lf seconds to thread \n",rank,temp);
            printf("The process took this time to finish threading %ld hours,%ld minutes,%ld seconds \n",t_hours,t_minutes,t_seconds);
            #endif
        }
        else{
            /*SCATTERING*/
            gettimeofday(&t_timer2_start, NULL); // start communication timer
            MPI_Scatterv(u_speed,sendcounts,displs,MPI_FLOAT,levels,GRID_POINTS*levels_per_proc,MPI_FLOAT,0,MPI_COMM_WORLD);
            gettimeofday(&t_timer2_finish, NULL);
            temp=time_diff(&t_timer2_start, &t_timer2_finish);
            #ifdef DEBUG
            convert_time_hour_sec(temp,&t_hours,&t_minutes,&t_seconds);
            printf("The processes %d took %lf seconds to scatter \n",rank,temp);
            printf("The process took this time to finish scattering %ld hours,%ld minutes,%ld seconds \n",t_hours,t_minutes,t_seconds);
            #endif
            /*THREADING*/
            gettimeofday(&t_timer2_start, NULL); // start communication timer
            threading(sum, levels, rank, size, levels_per_proc, thread_count);
            gettimeofday(&t_timer2_finish, NULL);
            temp=time_diff(&t_timer2_start, &t_timer2_finish);
            #ifdef DEBUG
            convert_time_hour_sec(temp,&t_hours,&t_minutes,&t_seconds);
            printf("The processes %d took %lf seconds to thread \n",rank,temp);
            printf("The process took this time to finish threading %ld hours,%ld minutes,%ld seconds \n",t_hours,t_minutes,t_seconds);
            #endif
        }

        gettimeofday(&t_timer2_start, NULL); // start communication timer
        MPI_Reduce(sum, final_averages,GRID_POINTS, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
        gettimeofday(&t_timer2_finish, NULL);
        temp=time_diff(&t_timer2_start, &t_timer2_finish);
        free(sum);
        gettimeofday(&t_timer1_finish, NULL);
        temp2=time_diff(&t_timer1_start, &t_timer1_finish);
        if(rank==0){
            #ifdef DEBUG
            convert_time_hour_sec(temp,&t_hours,&t_minutes,&t_seconds);
            printf("The processes %d took %lf seconds to reduce \n",rank,temp);
            printf("The process took to reduce %ld hours,%ld minutes,%ld seconds \n",t_hours,t_minutes,t_seconds);
            convert_time_hour_sec(temp2,&t_hours,&t_minutes,&t_seconds);
            printf("The processes %d took %lf seconds to finish 1 time data point \n",rank,temp2);
            printf("The process took this time to finish 1 time data point %ld hours,%ld minutes,%ld seconds \n",t_hours,t_minutes,t_seconds);
            #endif
            net_write(final_averages,k);
        }
        free(final_averages);
    }
    if (rank == 0){
        /*TIME END T3*/
        gettimeofday(&t_timer3_finish, NULL); 
        temp=time_diff(&t_timer3_start, &t_timer3_finish);
        convert_time_hour_sec(temp,&t_hours,&t_minutes,&t_seconds);
        printf("#########THE END OF COMPUTATION OVER ALL TIMESTEPS#######\n");
        printf("The time taken to paralize everything for all of the time steps %lf seconds\n",temp);
        printf("The time taken to paralize everything for all of the time steps %ld hours,%ld minutes,%ld seconds \n",t_hours,t_minutes,t_seconds);
    }
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
void threading(float * sum,float levels[][GRID_POINTS],int  rank,int size,int levels_per_proc,int  thread_count)
{
    int boundary;
    int i, j;
    if (((rank)-1) == (size))
    {
        boundary = (N_NZ1-(levels_per_proc * ((size)-1)));
    }
    else{
        boundary = (levels_per_proc);
    }
    for (i = 0; i < boundary;i++){
        for (j = 0; j < GRID_POINTS;j++){
                sum[j] += levels[i][j]/N_NZ1;
            }
    }
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