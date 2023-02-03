#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <netcdf.h>
#include <mpi.h>
#include <omp.h>
#include <sys/time.h>


/*MACROS START*/
#define FILE_NAME "writing.nc"
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); return 2;} //just simple macro to handle error with NETCDF
#define ERR_spec(e) {printf("Error: %s\n", nc_strerror(e));} //just simple macro to handle error with NETCDF
#define HEIGHT_1_NAME "height"  // the variable is called nz in the file. I knew it using ncdump command
#define TIME "time"  // the variable is called time in the file. I knew it using ncdump command
#define UNOD "speed"  // the variable is called unod in the file. I knew it using ncdump command
#define N_NZ1 10            // the variable nz1 has 69 entries
#define N_TIME 4            // the variable time has 12 entries
#define NDIMS 3
#define GRID_POINTS 5
// #define DEBUG 1
// for write
#define FILE_NAME2 "map_summarized_fuldataset_serial.nc"
#define UNITS_speed "m_s"
#define UNITS "units"
#define UNITS_time "s"
#define NDIMS_wr 3
#define SPLIT_COMM 2
/*MACROS end*/

double time_diff(struct timeval *start, struct timeval *end);
void convert_time_hour_sec( double seconds, long int *h, long int *m, long int *s);
void net_write(float * final_averages,int k);
void threading(float *sum, float **levels, int levels_per_proc);
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
    size_t start[NDIMS], count[NDIMS];

    /*Time variables to be used to see how much time each process takes*/

    struct timeval t_timer2_start;
    struct timeval t_timer2_finish;
    struct timeval t_timer3_start;
    struct timeval t_timer3_finish;
    double t_start = 0;
    double t_finish = 0;
    double t_threading = 0;
    double t_threading_sum;
    double t_reducing = 0;
    // double t_reducing_sum;
    long int t_seconds = 0;
    long int t_minutes = 0;
    long int t_hours = 0;
    double temp;
    /*DYNAMIC ALLOCATION*/
    int time_per_proc =ceil((double)N_TIME / SPLIT_COMM);//3
    /*color between 0 to 4*/
    int limit_time = (color + 1) * time_per_proc;

    if (limit_time > N_TIME)
    {
        limit_time = N_TIME;
    }

    int levels_per_proc = ceil((double)N_NZ1 / row_size);/*if 2.3 is passed to ceil(), it will return 3*/
    if ((retval = nc_open(FILE_NAME, NC_NOWRITE, &ncid)))
    ERR(retval);
    if ((retval = nc_inq_varid(ncid, UNOD, &unod_id)))
    ERR(retval);
    int count_levels_per_proc=levels_per_proc;
    if(row_rank==row_size-1){
        count_levels_per_proc = N_NZ1-levels_per_proc*row_rank;
        printf("%dfor the row rank equal to %d \n", count_levels_per_proc,row_rank);
    }
    count[0] = 1;/*1 time*/
    count[1] = count_levels_per_proc;/*14 level*/
    count[2] = GRID_POINTS;/*all gridpoints*/
    start[1] = row_rank*levels_per_proc;/*start from level 0 to 18 then 18-29*/
    start[2] = 0;
    printf("For row %d i have the following values %zu\n ",row_rank,start[1]);
    // float  u_speed[count_levels_per_proc][GRID_POINTS];
    // float **u_speed = (float**)malloc(sizeof(float *) * levels_per_proc);
    // for (i = 0; i < levels_per_proc; i++){
    //     u_speed[i] = (float *) calloc(GRID_POINTS ,sizeof(float));
    //     if (u_speed[i] == NULL) {
    //         printf("A Problem will occur now in allocate the u_speed");
    //     }
    // }
    // static float u_speed[18][GRID_POINTS];
    float **u_speed;
    float *p = calloc(count_levels_per_proc*GRID_POINTS,sizeof(float));
    (u_speed) = malloc(count_levels_per_proc*sizeof(float*));
    for (i=0; i<count_levels_per_proc; i++){
       (u_speed)[i] = &(p[i*GRID_POINTS]);
    }
    // printf("count level per process %d for process %d\n",color * time_per_proc,rank );
    /*START*/
    if (row_rank == 0){
            printf("Number of processes: %d (levels being read for each process: %d)\n", size, count_levels_per_proc);
    /*TIME START T3*/
            printf("#########THE START OF COMPUTATION OVER ALL TIMESTEPS#######\n");
            gettimeofday(&t_timer3_start, NULL); //start timer of rank0
    }
    
    for (k = color * time_per_proc;k< limit_time; k++){

        /*Instlalizat variables*/

        float * sum_u_speed;
        sum_u_speed = (float *)calloc(GRID_POINTS, sizeof(float));
        float *final_averages;
        final_averages =(float *)calloc(GRID_POINTS, sizeof(float));
        if (sum_u_speed == NULL || final_averages == NULL) {
            printf("A Problem will occur now ");
        }

        start[0] = k;
        t_threading_sum = 0;
        if ((retval = nc_get_vara_float(ncid, unod_id, start, count, &u_speed[0][0])))
            ERR(retval);
        // break;
        /*THREADING*/
        t_start = omp_get_wtime();
        // threading(sum_u_speed, u_speed, count_levels_per_proc);
        // if(row_rank==row_size-1){
        //     printf("the value is %d", count_levels_per_proc);
        // }
        for (i = 0; i < count_levels_per_proc; i++){
            for (j = 0; j < GRID_POINTS; j++){
                    sum_u_speed[j] += u_speed[i][j]/ N_NZ1;;
            }
        }
        t_finish = omp_get_wtime();
        t_threading = t_finish - t_start;
        #ifdef DEBUG
        convert_time_hour_sec(t_threading,&t_hours,&t_minutes,&t_seconds);
        printf("The processes %d took %lf seconds to thread \n",rank,t_threading);
        printf("The process took this time to finish threading %ld hours,%ld minutes,%ld seconds \n",t_hours,t_minutes,t_seconds);
        #endif



        /*REDUCE*/
        gettimeofday(&t_timer2_start, NULL); // start communication timer
        MPI_Reduce(sum_u_speed, final_averages,GRID_POINTS, MPI_FLOAT, MPI_SUM, 0, row_comm);
        MPI_Reduce(&t_threading, &t_threading_sum, 1, MPI_DOUBLE, MPI_SUM, 0, row_comm);
        gettimeofday(&t_timer2_finish, NULL);
        t_reducing=time_diff(&t_timer2_start, &t_timer2_finish);
        #ifdef DEBUG
        convert_time_hour_sec(t_reducing,&t_hours,&t_minutes,&t_seconds);
        printf("The processes %d took %lf seconds to reduce \n",rank,t_reducing);
        printf("The process took this time to finish reducing %ld hours,%ld minutes,%ld seconds \n",t_hours,t_minutes,t_seconds);
        #endif
        if(row_rank==0){
                    net_write_sep_files(final_averages,k);
                    // net_write(final_averages,k);
        }
        free(sum_u_speed);
        free(final_averages);
    }
    if (rank == 0){
        gettimeofday(&t_timer3_finish, NULL);
        /*TIME END 3*/
        temp=time_diff(&t_timer3_start, &t_timer3_finish);
        convert_time_hour_sec(temp,&t_hours,&t_minutes,&t_seconds);
        printf("The time taken to paralize everything for all of the files %lf seconds\n",temp);
        printf("The time taken to paralize everything for all of the files %ld hours,%ld minutes,%ld seconds \n",t_hours,t_minutes,t_seconds);
    }

    /*CLOSING FILE*/
    // for (i = 0; i < count_levels_per_proc; i++){
    //         free(u_speed[i]);      
    // }
    if ((retval = nc_close(ncid)))
        ERR(retval);
    // for (int i = 0; i < levels_per_proc; i++){
    //     free(u_speed[i]);
    // }
    // printf("This process has about to close right:on this node %s,with this ranking %d out of %d processes \n",
    //        processor_name, rank, size);
    MPI_Comm_free(&row_comm);
    MPI_Finalize();
    return 0;
}


void threading(float * sum,float ** levels,int levels_per_proc)
{
    int i, j;
    for (i = 0; i < levels_per_proc; i++){
        for (j = 0; j < GRID_POINTS; j++){
                    sum[j] += levels[i][j]/ N_NZ1;;
        }
    }

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