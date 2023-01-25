#include <omp.h>
#include<stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <mpi.h>

#define GRIDPOITN 8
#define level 4
/*
 gcc -g -Wall -fopenmp -o testing_idea testing_idea.c -lm
./testing_idea 5
*/
// echo |cpp -fopenmp -dM |grep -i open
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
void threading(int sum [level], int arr_dim[level][GRIDPOITN],int *thread_count );
void threading2(int sum[level], int arr_dim[level][GRIDPOITN], int *thread_count);
void threading3(int sum[level], int arr_dim[level][GRIDPOITN], int *thread_count);
void threading4(int sum[level], int arr_dim[level][GRIDPOITN], int *thread_count);
void threading5(int sum[level], int arr_dim[level][GRIDPOITN], int *thread_count);
int main(int argc, char *argv[])
{
    // int arr[8]  = {1, 2, 3, 4, 5, 6, 7, 8};
    int arr_dim[level][GRIDPOITN]  = {{1, 2, 3, 4, 5, 6, 7, 8},{1, 2, 3, 4, 5, 6, 7, 8},{1, 2, 3, 4, 5, 6, 7, 8},{1, 2, 3, 4, 5, 6, 7, 8}};
    // struct timeval t_timer1_start;/*timer for process 0*/
    // struct timeval t_timer1_finish;
    // long int t_hours, t_minutes, t_seconds;
    // double t_nc_reading_time = 0;
    int sum[GRIDPOITN]={0};
    int i;
    int threadnumb;
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
// int y = 0;
// int thread_count;
// thread_count = strtol(argv[1], NULL, 10);
#pragma omp parallel 
threadnumb = omp_get_num_threads();
    printf("The number of threads are %d \n",threadnumb);
    // gettimeofday(&t_timer1_start, NULL);
    threading2(sum,arr_dim,&threadnumb);
    threading3(sum,arr_dim,&threadnumb);
    threading4(sum,arr_dim,&threadnumb);
    threading5(sum,arr_dim,&threadnumb);
    // gettimeofday(&t_timer1_finish, NULL);
    // t_nc_reading_time = time_diff(&t_timer1_start, &t_timer1_finish);
    // convert_time_hour_sec(t_nc_reading_time,&t_hours,&t_minutes,&t_seconds);
    // printf("The time taken to do Nc read is %lf seconds\n",t_nc_reading_time);
    // printf("The time taken to do Nc read is %ld hours,%ld minutes,%ld seconds \n", t_hours,t_minutes,t_seconds);
    // printf("The output\n");
    for (i = 0; i < GRIDPOITN; i++){
        printf("%d,",sum[i]);
    }
    
    printf("\n");
     return 0;
}
void threading(int sum [level],int arr_dim [level][GRIDPOITN], int *thread_count){
    int i, j;
    double start, finish, elapsed;
    start = omp_get_wtime();
#pragma omp parallel for default(none) private(i, j)\
        shared(sum,arr_dim) 
    for (i = 0; i < level; i++){
        for (j = 0; j < GRIDPOITN; j++){
            sum[j] += arr_dim[i][j];
        }
    }
    finish = omp_get_wtime();
    elapsed = finish - start;
    printf("Elapsed time = %e seconds\n", elapsed);
}
void threading4(int sum [level],int arr_dim [level][GRIDPOITN], int *thread_count){
    int i, j;
    double start, finish, elapsed;
    start = omp_get_wtime();
            for (i = 0; i < level; i++){
                for (j = 0; j < GRIDPOITN; j++){
                    sum[j] += arr_dim[i][j];
                }
            }
    finish = omp_get_wtime();
    elapsed = finish - start;
    printf("Elapsed time for threading  without threds 2  = %e seconds\n", elapsed);

}
void threading5(int sum [level],int arr_dim [level][GRIDPOITN], int *thread_count){
    int i, j;
    double start, finish, elapsed;
    start = omp_get_wtime();
    
            for ( i = 0; i < GRIDPOITN; i++){
                for ( j = 0; j < level; j++) {
                    sum[i] += arr_dim[j][i];
                }
            }
    
    finish = omp_get_wtime();
    elapsed = finish - start;
    printf("Elapsed time for threading without threds 3  = %e seconds\n", elapsed);

}
void threading2(int sum [level],int arr_dim [level][GRIDPOITN], int *thread_count){
    int i, j;
    double start, finish, elapsed;
    start = omp_get_wtime();
    #pragma omp parallel
    {
        float S_private[GRIDPOITN] = {0};
        #pragma omp for  private(i, j) 
            for (i = 0; i < level; i++){
                for (j = 0; j < GRIDPOITN; j++){
                    S_private[j] += arr_dim[i][j];
                }
            }
        #pragma omp critical
        {
            for(int n=0; n<GRIDPOITN; ++n) {
            sum[n] += S_private[n];
            }
        }
    }
    finish = omp_get_wtime();
    elapsed = finish - start;
    printf("Elapsed time for threading 2  = %e seconds\n", elapsed);

}
void threading3(int sum [level],int arr_dim [level][GRIDPOITN], int *thread_count){
    int i, j;
    double start, finish, elapsed;
    start = omp_get_wtime();
    #pragma omp parallel
    {
        float S_private[GRIDPOITN] = {0};
        #pragma omp for  private(i, j) 
            for ( i = 0; i < GRIDPOITN; i++){
                for ( j = 0; j < level; j++) {
                    S_private[i] += arr_dim[j][i];
                }
            }
        #pragma omp critical
        {
            for(int n=0; n<GRIDPOITN; ++n) {
            sum[n] += S_private[n];
            }
        }
    }
    finish = omp_get_wtime();
    elapsed = finish - start;
    printf("Elapsed time for threading 3  = %e seconds\n", elapsed);

}