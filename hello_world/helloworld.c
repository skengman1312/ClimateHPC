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
 int arr[4][3]  = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
 // Print off a hello world message
 //printf("Hello world from process rank %d out of %d processors\n", world_rank, world_size);
 //   if (world_rank == 0)

 int * rec[3];
 int sendcnt = 3; /* how many items are sent to each process */
 int recvcnt = 3; /* how many items are received by each process */
 MPI_Scatter(arr, sendcnt, MPI_INT,
                &rec, recvcnt, MPI_INT, 0, MPI_COMM_WORLD);

 printf("My rank is %i and i got %i\n", world_rank, rec);


 //printf("%i\n", arr[1][world_rank]);

 // Finalize the MPI environment. No more MPI calls can be made after this
 MPI_Finalize();
 return 0;
}
