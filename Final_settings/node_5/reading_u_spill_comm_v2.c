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
#define FILE_NAME2 "map_summarized_fuldataset_complex.nc"
#define UNITS_speed "m_s"
#define UNITS "units"
#define UNITS_time "s"
#define NDIMS_wr 3
#define SPLIT_COMM 5
/*MACROS end*/

double time_diff(struct timeval *start, struct timeval *end);
void convert_time_hour_sec( double seconds, long int *h, long int *m, long int *s);
void net_write(float * final_averages,int k);
// void threading(float *sum, float **levels, int levels_per_proc);
void net_write_sep_files(float *final_averages, int k);

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
    /*SPILITING_COMMUNICATION*/
     // COMUNICATION settings    
    int color = rank / SPLIT_COMM; // Determine color based on row
    // original rank for ordering
    MPI_Comm row_comm;
    MPI_Comm_split(MPI_COMM_WORLD, color, rank, &row_comm);
    int row_rank, row_size;
    MPI_Comm_rank(row_comm, &row_rank);
    MPI_Comm_size(row_comm, &row_size);
    printf("greetings:  %s, row rank %d for color %d the orginal rank is%d out of %d processes\n",
           processor_name, row_rank,color,rank, size);
    /* VARIABLES DEFINE START*/
    /*LOOOPING variables*/
    // int rec;
    int i, j;
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
    struct timeval t_timer2_start;
    struct timeval t_timer2_finish;
    struct timeval starttime;
    struct timeval endtime;
    double t_start = 0;
    double t_finish = 0;
    double t_threading = 0;
    double t_threading_sum;
    double t_reducing = 0;
    long int t_seconds = 0;
    long int t_minutes = 0;
    long int t_hours = 0;
    double walltimes_start[2];
    double walltimes_end[2];
    double passing_time = 0;
    double nc_reading ; //variable containing the local sum of all elapsed time for reading
    double threading_time = 0;  //variable containing the local sum of all elapsed time for writing the sum matrix
    double sum_nc_reading;        // variable containing the output of reduce operation, collecting all reading times of different processes
    double sum_threading_time; // variable containing the output of reduce operation, collecting all elaboration times of different processes
    double total_time_reading=0;
    double total_time_threading=0;
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

    int time_per_proc =ceil((double)N_TIME / SPLIT_COMM);//3
    /*color between 0 to 4*/
    int limit_time = (color + 1) * time_per_proc;

    if (limit_time > N_TIME)
    {
        limit_time = N_TIME;
    }

    if ((retval = nc_open(FILE_NAME, NC_NOWRITE, &ncid)))
    ERR(retval);
    if ((retval = nc_inq_varid(ncid, UNOD, &unod_id)))
    ERR(retval);
    float * u_speed;
    u_speed = (float *)calloc(GRID_POINTS, sizeof(float));
    int levels_per_proc = ceil((double)N_NZ1 / row_size);/*if 2.3 is passed to ceil(), it will return 3*/
    int limit = (row_rank + 1) * levels_per_proc;
    if (limit > N_NZ1)
    {
        limit = N_NZ1;
    }
    if (row_rank == 0){
        walltimes_start[1] = MPI_Wtime();
    }
    for (k = color * time_per_proc;k< limit_time; k++){
        /*Instlalizat variables*/
        float * sum_u_speed;
        sum_u_speed = (float *) calloc(GRID_POINTS, sizeof(float));
        float *final_averages;
        final_averages =(float *) calloc(GRID_POINTS, sizeof(float));
        if (sum_u_speed == NULL || final_averages == NULL) {
            printf("A Problem will occur now ");
        }       
        start[0] = k;
        t_threading_sum = 0;
        nc_reading = 0;
        threading_time = 0;

        /*THREADING*/
        t_start = omp_get_wtime();
        count[0] = 1;/*1 time*/
        count[1] = 1;/*1 level*/
        count[2] = GRID_POINTS;/*all gridpoints*/
        start[2] = 0;
        for (j = row_rank * levels_per_proc;j< limit; j++){
            start[1] = j;
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
        t_finish = omp_get_wtime();
        t_threading = t_finish - t_start;
        MPI_Reduce(&nc_reading, &sum_nc_reading, 1, MPI_DOUBLE, MPI_SUM, 0, row_comm);
        MPI_Reduce(&threading_time, &sum_threading_time, 1, MPI_DOUBLE, MPI_SUM, 0, row_comm);

        /*REDUCE*/
        gettimeofday(&t_timer2_start, NULL); // start communication timer
        MPI_Reduce(sum_u_speed, final_averages,GRID_POINTS, MPI_FLOAT, MPI_SUM, 0, row_comm);
        MPI_Reduce(&t_threading, &t_threading_sum, 1, MPI_DOUBLE, MPI_SUM, 0, row_comm);
        gettimeofday(&t_timer2_finish, NULL);
        t_reducing=time_diff(&t_timer2_start, &t_timer2_finish);

        if(row_rank==0){
                // net_write_sep_files(final_averages,k);
                total_time_reading += sum_nc_reading / row_size;
                total_time_threading += sum_threading_time / row_size;
                walltimes_end[0] = MPI_Wtime();  
                convert_time_hour_sec(t_reducing,&t_hours,&t_minutes,&t_seconds);
                printf("##### THE BEGINING OF THE RESULT OF INSTANCE %d ##### \n",k);
                printf("The time taken to do 3 MPI REDUCE is %lfseconds \n",t_reducing);
                printf("The time taken to do 3 MPI REDUCE is %ld hours,%ld minutes,%ld seconds \n",t_hours,t_minutes,t_seconds);
                printf("Average reading time per process: %.7f\n", sum_nc_reading / row_size);
                printf("Average threading time per process: %.7f\n\n", sum_threading_time / row_size);
                printf("TotaL looping time for %lf seconds for a number rows per column equal to %d \n",t_threading_sum/row_size,SPLIT_COMM);
                printf("##### THE END OF THE RESULT OF INSTANCE %d #####\n ",k);
                printf("\n");
                // net_write(final_averages,k);
        }
        free(sum_u_speed);
        free(final_averages);
    }
    if (row_rank == 0){
        walltimes_end[1] = MPI_Wtime();
        /*TIME END 3*/
         printf("####THIS IS THE FINAL RESULTSSSS###\n");
        printf("Average reading time per process over X month: %.7f\n", total_time_reading / time_per_proc);
        printf("Average threading time per process over X month: %.7f\n\n", total_time_threading / time_per_proc);
        printf("Total walltime from the begining till the end for ALL time step is %lf seconds \n", walltimes_end[1] - walltimes_start[1]);
    }

    /*CLOSING FILE*/
    // /* free the pointers into the memory */
    free(u_speed);
    if ((retval = nc_close(ncid)))
        ERR(retval);
    MPI_Comm_free(&row_comm);
    MPI_Finalize();
    return 0;
}

void net_write_sep_files(float * final_averages, int k){
    /*START creating file */
    int retval, ncid2;
    int time_dimid; // time dimension id 
    int speed_dimid; // speed dimension id 
    int dimid[2]; // for wrtieing 
    int var_speed_id;// variable speed id
    char *pointer;
    switch (k)
    {
        case 0:
            pointer="nc_file_January.nc";
            break;
        case 1: 
            pointer="nc_file_February.nc";
            break;
        case 2: 
            pointer="nc_file_March.nc";
            break;
        case 3: 
            pointer="nc_file_April.nc";
            break;
        case 4: 
            pointer="nc_file_May.nc";
            break;
        case 5: 
            pointer="nc_file_June.nc";
            break;
        case 6: 
            pointer="nc_file_July.nc";
            break;
        case 7: 
            pointer="nc_file_August.nc";
            break;
        case 8: 
            pointer="nc_file_September.nc";
            break;
        case 9: 
            pointer="nc_file_October.nc";
            break;
        case 10: 
            pointer="nc_file_November.nc";
            break;
        default: 
            pointer="nc_file_December.nc";
    }
    if ((retval = nc_create(pointer, NC_CLOBBER, &ncid2))) // ncclober to overwrite the file
        ERR_spec(retval);
    if ((retval = nc_def_dim(ncid2, TIME, NC_UNLIMITED, &time_dimid)))
            ERR_spec(retval);
    if ((retval = nc_def_dim(ncid2, UNOD, GRID_POINTS, &speed_dimid)))
            ERR_spec(retval);
    dimid[0]=time_dimid;
    dimid[1]=speed_dimid;
    if ((retval = nc_def_var(ncid2, "speed", NC_FLOAT, 2, dimid, &var_speed_id)))// define the varibael
            ERR_spec(retval);
    if ((retval = nc_put_att_text(ncid2, var_speed_id, UNITS, 
				 strlen(UNITS_speed), UNITS_speed)))
            ERR_spec(retval);
        /* End define mode. */
    if ((retval = nc_enddef(ncid2)))
            ERR_spec(retval);
    size_t start_1[2]={0,0};
    size_t count_1[2]={1,GRID_POINTS};
    if ((retval = nc_put_vara_float(ncid2,var_speed_id, start_1, count_1,&final_averages[0])))ERR_spec(retval);
    if ((retval = nc_close(ncid2)))ERR_spec(retval);
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