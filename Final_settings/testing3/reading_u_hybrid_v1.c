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
#define FILE_NAME2 "map_summarized_fuldataset_hybrid_v1.nc"
#define UNITS_speed "m_s"
#define UNITS "units"
#define UNITS_time "s"
#define NDIMS_wr 3
// #define SPLIT_COMM 4
/*MACROS end*/

double time_diff(struct timeval *start, struct timeval *end);
void convert_time_hour_sec( double seconds, long int *h, long int *m, long int *s);
void net_write(float * final_averages,int k);
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
    printf("greetings:  %s, rank %d out of %d processes\n",processor_name, rank, size);
    int rec;
    int i,k;
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
    
    struct timeval t_timer1_start;/*timer for process 0*/
    struct timeval t_timer1_finish;
    struct timeval t_timer2_start;
    struct timeval t_timer2_finish;
    struct timeval t_timer3_start;
    struct timeval t_timer3_finish;
    double t_threading_reading_time ;
    double t_threading_reading_time_sum ;
    double t_threading_reading_time_Totalsum ;
    double t_comm_time ;
    double t_time_from_start;
    long int t_seconds = 0;
    long int t_minutes = 0;
    long int t_hours = 0;
    double temp;
 
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
    // start[1] = 0;
    start[2] = 0;

    /* START timer 3 to have calculate total time*/
    if (rank == 0){
        gettimeofday(&t_timer3_start, NULL); //start timer of rank0
    }
    for (k = 0;k<N_TIME; k++){
        start[0] = k;
        t_threading_reading_time = 0 ;
        t_threading_reading_time_sum = 0 ;
        t_threading_reading_time_Totalsum = 0;
        t_comm_time = 0 ;
        t_time_from_start = 0;
        float * sum_u_speed;
        sum_u_speed = (float*)calloc(GRID_POINTS, sizeof(float));
        float *final_averages;
        final_averages = (float *)calloc(GRID_POINTS, sizeof(float));
        if (rank == 0){
            printf("Number of processes: %d (levels being read for each process: %d)\n", size, levels_per_proc);
            /*TIME START T1*/
            gettimeofday(&t_timer1_start, NULL); //start timer of rank0
        }
        for (rec = rank * levels_per_proc; rec < limit; rec++){

            /*TIME START T2*/
            start[1] = rec;
            if ((retval = nc_get_vara_float(ncid, unod_id, start, count, &u_speed[0])))
                ERR(retval);



            /*TIME START T2*/
            gettimeofday(&t_timer2_start, NULL); //start reading timer
            #pragma omp for  private(i) schedule(guided)
            for (i = 0; i < GRID_POINTS;i++){
                sum_u_speed[i] += u_speed[i]/ N_NZ1;
            }
            /*
            schedule(static)
            schedule(dynamic)
            schedule(guided)
            */
            gettimeofday(&t_timer2_finish, NULL);
            /*TIME END T2*/
            
            /* calculate the time take for SUMMING*/
            t_threading_reading_time = time_diff(&t_timer2_start, &t_timer2_finish);
            t_threading_reading_time_sum+=t_threading_reading_time;
        }

        /*TIME START T2*/
        gettimeofday(&t_timer2_start, NULL); // start communication timer
        MPI_Reduce(&t_threading_reading_time_sum, &t_threading_reading_time_Totalsum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(sum_u_speed, final_averages,GRID_POINTS, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);    
        gettimeofday(&t_timer2_finish, NULL);
        /*TIME END T2*/
        t_comm_time=time_diff(&t_timer2_start, &t_timer2_finish);
        if (rank == 0){ 
            gettimeofday(&t_timer1_finish, NULL); //start timer of rank0
            t_time_from_start=time_diff(&t_timer1_start, &t_timer1_finish);
            /*TIME END T1*/
            net_write(final_averages,k);
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
        free(sum_u_speed);
        free(final_averages);
    }
    if (rank == 0){
        gettimeofday(&t_timer3_finish, NULL);
        /*TIME END 3*/
        temp=time_diff(&t_timer3_start, &t_timer3_finish);
        convert_time_hour_sec(temp,&t_hours,&t_minutes,&t_seconds);
        printf("####THIS IS THE FINAL RESULTSSSS###\n");
        printf("The time taken to paralize everything for all of the files %lf seconds\n",temp);
        printf("The time taken to paralize everything for all of the files %ld hours,%ld minutes,%ld seconds \n",t_hours,t_minutes,t_seconds);
    }
    free(u_speed);
    if ((retval = nc_close(ncid)))
        ERR(retval);
    printf("This process has about to close right:on this node %s,with this ranking %d out of %d processes \n",
           processor_name, rank, size);
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
    // if ((retval = nc_close(ncid2)))
    //         ERR(retval);
    size_t start_1[2]={0,0};
    size_t count_1[2]={1,GRID_POINTS};
    // if ((retval = nc_open(FILE_NAME2, NC_WRITE, &ncid)))ERR_spec(retval);
    // if ((retval = nc_inq_varid(ncid2, "speed", &unod_id)))ERR_spec(retval);
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