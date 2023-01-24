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
#define SUBARRAY_SIZE 2

// #define DEBUG 0
// #define ALL_procc 0

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
void threading(float * sum, float ** levels, int rank, int size, int levels_per_proc);
void net_write(float * final_averages,int k);
/*FUNCTIONS end*/
int main (int argc, char *argv[]){
    /* MPI  inizialization */
    MPI_Init(&argc, &argv);

    // MPI_Init(&argc, &argv);

    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
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
    double temp;
    double temp2;
    // static float u_speed[N_NZ1][GRID_POINTS] = {{0}};
    float **u_speed = malloc(sizeof(float *) * N_NZ1);
    for (i = 0; i < N_NZ1; i++)
        u_speed[i] = (float *)calloc(GRID_POINTS ,sizeof(float));
    /*Time variables to be used to see how much time each process takes*/
    struct timeval t_timer1_start;/*timer for process 0*/
    struct timeval t_timer1_finish;
    struct timeval t_timer2_start;
    struct timeval t_timer2_finish;
    struct timeval t_timer3_start;
    struct timeval t_timer3_finish;
    double t_start = 0;
    double t_finish = 0;
    double t_nc_scatter = 0;
    double t_threading= 0;
    double t_nc_scatter_sum;
    double t_threading_sum;
    long int t_seconds = 0;
    long int t_minutes = 0;
    long int t_hours = 0;
    // FILE *fpt;

    // CREATE WRITE FILE 
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
    /*START*/


     // COMUNICATION settings    
    int color = rank / 5; // Determine color based on row
    // original rank for ordering
    MPI_Comm row_comm;
    MPI_Comm_split(MPI_COMM_WORLD, color, rank, &row_comm);
    int row_rank, row_size;
    MPI_Comm_rank(row_comm, &row_rank);
    MPI_Comm_size(row_comm, &row_size);
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    // Print off a hello world message with processor name
    // printf("greetings:  %s, rank %d out of %d processes\n",
    //        processor_name, rank, size);
    printf("WORLD RANK/SIZE: %d/%d \t ROW RANK/SIZE: %d/%d\n",rank, size, row_rank, row_size);
    
    int time_per_proc =ceil((double)N_TIME / 5);
    int limit = (row_rank + 1) * time_per_proc;

    if (limit > N_TIME)
    {
        limit = N_TIME;
    }

    /* Opening file in all the process*/
    if(0==row_rank){

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
    int levels_per_proc = ceil((double)N_NZ1 / row_size);/*if 2.3 is passed to ceil(), it will return 3*/
    int sendcounts[row_size];// for scatterv
    int displs[row_size];// for scatterv
    int sizes[2] = {N_NZ1, GRID_POINTS};
    int subsizes[2] = {levels_per_proc, GRID_POINTS};
    int starts[2] = {0, 0};
    MPI_Datatype subarray,resizedtype;
    MPI_Type_create_subarray(SUBARRAY_SIZE,sizes,subsizes,starts,MPI_ORDER_C,MPI_FLOAT,&subarray);
    MPI_Type_create_resized(subarray, 0, GRID_POINTS*levels_per_proc*sizeof(float), &resizedtype);
    MPI_Type_commit(&resizedtype);
    if (row_rank == 0){
            printf("Number of processes: %d (levels being read for each process: %d)\n", row_size, levels_per_proc);
/*TIME START T3*/
            printf("#########THE START OF COMPUTATION OVER ALL TIMESTEPS#######\n");
            gettimeofday(&t_timer3_start, NULL); //start timer of rank0
    }
    // intializing the arrays
    for (i = 0; i < row_size;i++){
        sendcounts[i]= levels_per_proc * GRID_POINTS;
        displs[i] = i * levels_per_proc * GRID_POINTS;
    }
    sendcounts[row_size-1] = (N_NZ1-(levels_per_proc*(row_size-1)))*GRID_POINTS;
    // float levels[levels_per_proc][GRID_POINTS];
    float **levels = malloc(sizeof(float *) * levels_per_proc);
    for (i = 0; i < levels_per_proc; i++){
        levels[i] = (float *)calloc(GRID_POINTS ,sizeof(float));
    }
    for (k = row_rank * time_per_proc;k< limit; k++){
            float * sum;
            float *final_averages;
            sum = (float *)calloc(GRID_POINTS, sizeof(float));
            if (!sum)
            {
                fprintf(stderr, "Could not allocate summm %d\n", row_rank); // I am trying to send sthg that is nully so it will raise rror
                MPI_Abort(row_comm, 1);
            }
            final_averages = (float *)calloc(GRID_POINTS, sizeof(float));
            if (!final_averages)
            {
                fprintf(stderr, "Could not allocate finale averages %d\n", row_rank); // I am trying to send sthg that is nully so it will raise rror
                MPI_Abort(row_comm, 1);
            }
            t_nc_scatter_sum= 0;
            t_threading_sum= 0;

            if (row_rank == 0)
            {
                gettimeofday(&t_timer1_start, NULL); //start timer of rank0

            }
            if(row_rank==0){
                start[0] = k;
                if ((retval = nc_get_vara_float(ncid, unod_id, start, count, &u_speed[0][0])))
                            ERR(retval);
                printf("I read %lf, from process %d \n",u_speed[0][150],row_rank);
            }
            /*SCATTERING*/
            // MPI_Barrier(row_comm);
            gettimeofday(&t_timer2_start, NULL); // start communication timer
            MPI_Scatterv(&(u_speed[0][0]),sendcounts,displs,resizedtype,&(levels[0][0]),GRID_POINTS*levels_per_proc,MPI_FLOAT,0,row_comm);
            gettimeofday(&t_timer2_finish, NULL);
            t_nc_scatter=time_diff(&t_timer2_start, &t_timer2_finish);
            printf("i read %lf in row %d\n",levels[0][0],row_rank);
            free(final_averages);
            free(sum);
            break;
    }
    // fclose(fpt);
    for (i = 0; i < levels_per_proc; i++){
        free(levels[i]);
    }
    free(levels);
    for (i = 0; i < N_NZ1; i++){
        free(u_speed[i]);
    }
    free(u_speed);
    MPI_Comm_free(&row_comm);
    MPI_Type_free(&resizedtype);
    printf("Terminated properly\n");
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
    printf("%d,%d,%d",rank,size,levels_per_proc);
    // start = omp_get_wtime();

    for (i = 0; i < boundary; i++)
    {
        for (j = 0; j < GRID_POINTS; j++)
        {
            sum[j] += levels[i][j];
        }
            }
    // finish = omp_get_wtime();
    // elapsed = finish - start;
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