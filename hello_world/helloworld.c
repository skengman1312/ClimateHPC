#include<mpi.h>
#include<stdio.h>
/*
 mpicc -g -Wall -o mpi_peitro.out helloworld.c
 chmod u+x mpi_peitro 
 qsub test_helloworld.sh
*/
int main(int argc,char**argv){
 // Initialize the MPI environment. The two arguments to MPI Init are not
 // currently used by MPI implementations, but are there in case future
 // implementations might need the arguments.
 MPI_Init(NULL,NULL);

 // get the number of processes
 int world_size;
 MPI_Comm_size(MPI_COMM_WORLD,&world_size);

 // get the rank of processes
 int world_rank;
 MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
 int arr[8][3]  = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
                   13, 14, 15, 16, 17 ,18, 19 ,20, 21, 22, 23, 24};
 // Print off a hello world message
 //printf("Hello world from process rank %d out of %d processors\n", world_rank, world_size);
 //   if (world_rank == 0)

 int rec[2][3] = {0};
 float loc_avg[3] = {0};
 float avg[3] ={0};
 int sendcnt = 6; /* how many items are sent to each process */
 int recvcnt = 6; /* how many items are received by each process */
 MPI_Scatter(arr, sendcnt, MPI_INT,
                rec, recvcnt, MPI_INT, 0, MPI_COMM_WORLD);

 for (int i = 0; i < 3; i++){
     for (int j = 0; j < 2; ++j) {
         loc_avg[i] += rec[j][i];
     }
 }
    for (int i = 0; i < 3; ++i) {
        loc_avg[i] = loc_avg[i] / 2;
    }


 printf("My rank is %i and I recived from  %i to %i\n", world_rank, rec[0][0],rec[1][2]);
 printf("My loc_avg array is %g, %g, %g\n", loc_avg[0], loc_avg[1], loc_avg[2]);


 MPI_Reduce(loc_avg, avg, 3, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
 //printf("%i\n", arr[1][world_rank]);
    for (int i = 0; i < 3; ++i) {
        avg[i] = avg[i]/world_size;
    }
 if (world_rank == 0)
     printf("I am proc 0 and the collected avg is %g, %g, %g\n", avg[0], avg[1], avg[2]);
// TODO: Implement relative dimensions and test on real data
 // Finalize the MPI environment. No more MPI calls can be made after this
 MPI_Finalize();
 return 0;
}
