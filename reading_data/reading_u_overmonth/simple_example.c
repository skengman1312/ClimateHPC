#include <omp.h>
#include<stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>


/*
mpicc -std=c99 -g -Wall -fopenmp  -o testing_idea_2.out testing_idea_2.c -lm 

 gcc -g -Wall -fopenmp -o testing_idea testing_idea.c -lm
./testing_idea 5
*/
// echo |cpp -fopenmp -dM |grep -i open
int main(int argc, char* argv[]){
    MPI_Init(&argc, &argv);
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int arr_dim[4][8]  = {{1, 2, 3, 4, 5, 6, 7, 8},{-1, 2, -3, 4, -5, 6, -7, 8},{1, 2, 3, 4, 5, 6, 7, 8},{1, 2, 3, 4, 5, 6, 7, 8}};
    int levels_per_proc = ceil((double)4 / size);/*if 2.3 is passed to ceil(), it will return 3*/
    int limit = (rank + 1) * levels_per_proc;
    int rec;
    int j;
    int i;
    int thread_count;
#pragma omp parallel
    thread_count = omp_get_num_threads();
    if (limit > 4)
    {
        limit = 4;
    }
    static int sum_u_speed[8] = {0};
    int final_averages[8] = {0};
    /*
    you have 2 processes
    it takes first 2 arrays  and then take the second 2 arrays
    in the end of loop i would have added the two arrays togather and i have 1 insgle array 
    */
    for (rec = rank * levels_per_proc; rec < limit; rec++){
        # pragma omp parallel for num_threads(thread_count) private(i) 
        for (i = 0; i < 8; i++){
            sum_u_speed[i] += arr_dim[rec][i];
            int my_rank = omp_get_thread_num();// get the number of the thread
            printf("Hello from thread %d containing %d from processes %d\n", my_rank, sum_u_speed[i],rank);
        }
    }
    MPI_Reduce(sum_u_speed, final_averages,8, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank==0){
        printf("It worked and this is the output\n");
        for (j = 0; j < 8; j++){
            printf("%d\t", final_averages[j]);
        }
    }
     return 0;
 }
