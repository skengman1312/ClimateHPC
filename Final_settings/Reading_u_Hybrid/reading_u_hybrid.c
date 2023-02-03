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
// #define DEBUG 1
// for write
#define FILE_NAME2 "map_summarized_fuldataset.nc"
#define UNITS_speed "m_s"
#define UNITS "units"
#define UNITS_time "s"
#define NDIMS_wr 3
/*MACROS end*/


double time_diff(struct timeval *start, struct timeval *end);
void convert_time_hour_sec( double seconds, long int *h, long int *m, long int *s);
void net_write(float * final_averages,int k);
int malloc2D(float ***array, int n, int m);
int free2D(float ***array);
void threading(float *sum, float **levels, int rank, int size, int levels_per_proc);
void Check_for_error(int local_ok, char fname[], char message[],
                     MPI_Comm comm);

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
    printf("greetings:  %s, rank %d out of %d processes\n",processor_name, rank, size);
    /* VARIABLES DEFINE START*/
    int ncid;
    int unod_id;
    int retval;
    /*Netcdf id for write*/
    int ncid2; // write file nmae
    int time_dimid; // time dimension id 
    int speed_dimid; // speed dimension id 
    int dimid[2];
    int var_speed_id;// variable speed id 
    size_t start[NDIMS], count[NDIMS];//used for reading


    /*LOOOPING variables*/
    // int rec;
    int i;
    int k;
    int levels_per_proc = ceil((double)N_NZ1 / size);/*if 2.3 is passed to ceil(), it will return 3*/
    int sendcounts[size];// for scatterv
    int displs[size];// for scatterv
    double temp;
    double temp2;
    // intializing the arrays
    for (i = 0; i < size;i++){
        sendcounts[i]= levels_per_proc * GRID_POINTS;
        displs[i] = i * levels_per_proc * GRID_POINTS;
    }
    sendcounts[size-1] = (N_NZ1-(levels_per_proc*(size-1)))*GRID_POINTS;
    static float u_speed[N_NZ1][GRID_POINTS] = {{0}};
    float levels[15][GRID_POINTS]={{0}} ;

    // float **u_speed;
    // malloc2D(&u_speed,N_NZ1,GRID_POINTS);
    // float **levels;
    // malloc2D(&levels,levels_per_proc,GRID_POINTS);
    // float **levels = (float**) malloc(sizeof(float *) * levels_per_proc);
    // float **levels;
    // for (int i = 0; i < levels_per_proc; i++){
    //     levels[i] = (float *)calloc(GRID_POINTS ,sizeof(float));
    // }
    // malloc2D(&levels,levels_per_proc,GRID_POINTS);
    // float levels[levels_per_proc][GRID_POINTS];
    /*Time variables to be used to see how much time each process takes*/
    struct timeval t_timer1_start;/*timer for process 0*/
    struct timeval t_timer1_finish;
    struct timeval t_timer2_start;
    struct timeval t_timer2_finish;
    struct timeval t_timer3_start;
    struct timeval t_timer3_finish;
    double t_start = 0;
    double t_finish = 0;
    double t_scatter = 0;
    double t_threading= 0;
    double t_scatter_sum;
    double t_threading_sum;
    long int t_seconds = 0;
    long int t_minutes = 0;
    long int t_hours = 0;

    
    /*Creating the file */
    
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
        /*Process 0 is scatter everything*/
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
    if (rank == 0){
            printf("Number of processes: %d (levels being read for each process: %d)\n", size, levels_per_proc);
    /*TIME START T3*/
            printf("#########THE START OF COMPUTATION OVER ALL TIMESTEPS#######\n");
            gettimeofday(&t_timer3_start, NULL); //start timer of rank0
    }



    for (k = 0;k<N_TIME; k++){
        float * sum;
        float *final_averages;
        sum = (float *)calloc(GRID_POINTS, sizeof(float));
        final_averages = (float *)calloc(GRID_POINTS, sizeof(float));
        t_scatter_sum= 0;
        t_threading_sum= 0;
        if (rank == 0)
        {
            gettimeofday(&t_timer1_start, NULL); 
        }
        /*I am reading 69*8m float array*/
        if(rank==0){
            start[0] = k;
            if ((retval = nc_get_vara_float(ncid, unod_id, start, count, &u_speed[0][0])))
                        ERR(retval);
        }
        gettimeofday(&t_timer2_start, NULL); // start communication timer
        if(rank==0){
            MPI_Scatterv(u_speed,sendcounts,displs,MPI_FLOAT,levels,GRID_POINTS*levels_per_proc,MPI_FLOAT,0,MPI_COMM_WORLD);
        }
        else{
            MPI_Scatterv(u_speed,sendcounts,displs,MPI_FLOAT,levels,GRID_POINTS*levels_per_proc,MPI_FLOAT,0,MPI_COMM_WORLD);
        }
        gettimeofday(&t_timer2_finish, NULL);
        t_scatter=time_diff(&t_timer2_start, &t_timer2_finish);
        #ifdef DEBUG
        convert_time_hour_sec(t_scatter,&t_hours,&t_minutes,&t_seconds);
        printf("The processes %d took %lf seconds to scatter \n",rank,t_scatter);
        printf("The process took this time to finish scattering %ld hours,%ld minutes,%ld seconds \n",t_hours,t_minutes,t_seconds);
        #endif

        break;

        // /*THREADING*/
        // t_start = omp_get_wtime();
        // threading(sum, levels, rank, size, levels_per_proc);
        // t_finish = omp_get_wtime();
        // t_threading = t_finish - t_start;
        // #ifdef DEBUG
        // convert_time_hour_sec(t_threading,&t_hours,&t_minutes,&t_seconds);
        // printf("The processes %d took %lf seconds to thread \n",rank,t_threading);
        // printf("The process took this time to finish threading %ld hours,%ld minutes,%ld seconds \n",t_hours,t_minutes,t_seconds);
        // #endif

        // gettimeofday(&t_timer2_start, NULL); // start communication timer
        // MPI_Reduce(sum, final_averages,GRID_POINTS, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
        // MPI_Reduce(&t_scatter, &t_scatter_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        // MPI_Reduce(&t_threading, &t_threading_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        // gettimeofday(&t_timer2_finish, NULL);
        // temp=time_diff(&t_timer2_start, &t_timer2_finish);
        // free(sum);
        // temp2=time_diff(&t_timer1_start, &t_timer1_finish);
        // if(rank==0){
        //     t_scatter_sum /= size;
        //     t_threading_sum /= size;
        //     convert_time_hour_sec(temp,&t_hours,&t_minutes,&t_seconds);
        //     printf("The processes %d took %lf seconds to reduce \n",rank,temp);
        //     printf("The process took to reduce %ld hours,%ld minutes,%ld seconds \n",t_hours,t_minutes,t_seconds);
        //     convert_time_hour_sec(temp2,&t_hours,&t_minutes,&t_seconds);
        //     printf("The processes %d took %lf seconds to finish 1 time data point \n",rank,temp2);
        //     printf("The process took this time to finish 1 time data point %ld hours,%ld minutes,%ld seconds \n",t_hours,t_minutes,t_seconds);
        //     net_write(final_averages,k);
        // }
        // free(final_averages);
    }

    if (rank == 0){
        /*TIME END T3*/
        gettimeofday(&t_timer3_finish, NULL); 
        temp=time_diff(&t_timer3_start, &t_timer3_finish);
        convert_time_hour_sec(temp,&t_hours,&t_minutes,&t_seconds);
        // fprintf(fpt,"0, 0, 0, 0, %lf\n",temp);        
        printf("#########THE END OF COMPUTATION OVER ALL TIMESTEPS#######\n");
        printf("The time taken to paralize everything for all of the time steps %lf seconds\n",temp);
        printf("The time taken to paralize everything for all of the time steps %ld hours,%ld minutes,%ld seconds \n",t_hours,t_minutes,t_seconds);
    }
    MPI_Finalize();
    return 0;
}
/*-------------------------------------------------------------------*/
void Check_for_error(
      int       local_ok   /* in */, 
      char      fname[]    /* in */,
      char      message[]  /* in */, 
      MPI_Comm  comm       /* in */) {
   int ok;

   MPI_Allreduce(&local_ok, &ok, 1, MPI_INT, MPI_MIN, comm);
   if (ok == 0) {
      int my_rank;
      MPI_Comm_rank(comm, &my_rank);
      if (my_rank == 0) {
         fprintf(stderr, "Proc %d > In %s, %s\n", my_rank, fname, 
               message);
         fflush(stderr);
      }
      MPI_Finalize();
      exit(-1);
   }
}  /* Check_for_error */




int free2D(float ***array) {
    /* free the memory - the first element of the array is at the start */
    free(&((*array)[0][0]));

    /* free the pointers into the memory */
    free(*array);

    return 0;
}
int malloc2D(float ***array, int n, int m) {
    int i,local_ok;
    /* allocate the n*m contiguous items */
    float *p = calloc(n*m,sizeof(float));

    if (!p) return -1;

    /* allocate the row pointers into the memory */
    (*array) = malloc(n*sizeof(float*));
    if (!(*array)) {
       free(p);
       return -1;
    }
    if ( ( p == NULL ) || ( array == NULL ) ) {
            local_ok = 0;
            
    }
    Check_for_error(local_ok, "Main for loop", "Can't allocate memory vector", MPI_COMM_WORLD);

    /* set up the pointers into the contiguous memory */
    for (i=0; i<n; i++)
       (*array)[i] = &(p[i*m]);

    return 0;
}
void threading(float * sum,float ** levels,int  rank,int size,int levels_per_proc)
{
    int boundary;
    int i, j;
    // double start, finish, elapsed;

    if (((rank)-1) == (size))
    {
        boundary = (N_NZ1-(levels_per_proc * ((size)-1)));
    }
    else{
        boundary = (levels_per_proc);
    }
    // start = omp_get_wtime();
   #pragma omp parallel
    {
        float *S_private;
        S_private = (float *)calloc(GRID_POINTS, sizeof(float));
        #pragma omp for  private(i, j) 
            for (i = 0; i < boundary; i++){
                for (j = 0; j < GRID_POINTS; j++){
                    S_private[j] += levels[i][j];
                }
            }
        #pragma omp critical
        {
            for(int n=0; n<GRID_POINTS; ++n) {
                sum[n] += S_private[n];
            }
        }
        free(S_private);
    }
    // finish = omp_get_wtime();
    // elapsed = finish - start;
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