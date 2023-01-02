#include<mpi.h>
#include<stdio.h>
/*
 mpicc -g -Wall -o mpi_peitro helloworld.c
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

 // Print off a hello world message
 printf("Hello world from process rank %d out of %d processors\n", world_rank, world_size);


 // Finalize the MPI environment. No more MPI calls can be made after this
 MPI_Finalize();
 return 0;
}
