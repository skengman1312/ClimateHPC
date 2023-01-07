#include <omp.h>
#include<stdio.h>
#include <stdlib.h>
/*
 gcc -g -Wall -fopenmp -o testing_idea testing_idea.c -lm
./testing_idea 5
*/
// echo |cpp -fopenmp -dM |grep -i open
int main(int argc, char* argv[]){
    int arr[8]  = {1, 2, 3, 4, 5, 6, 7, 8};
    int i;
    int sum[8];
    // int y = 0;
    int thread_count;
    thread_count = strtol(argv[1], NULL, 10);
# pragma omp parallel for num_threads(thread_count) private(i) 
    for (i = 0; i < 8; i++){
        sum[i] += arr[i];
        int my_rank = omp_get_thread_num();// get the number of the thread
        printf("Hello from thread %d containing %d\n", my_rank, sum[i]);
    }

     return 0;
 }
