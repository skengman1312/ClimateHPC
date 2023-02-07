#include <stdio.h>
#include <string.h>
#include <netcdf.h>
#include<mpi.h>
#include <sys/time.h>

int main () {
    // Initialize the MPI environment. The two arguments to MPI Init are not
    // currently used by MPI implementations, but are there in case future
    // implementations might need the arguments.
    MPI_Init(NULL, NULL);

    // /* VARIABLES DEFINE START*/

    // Get the rank and size in the original communicator
    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int color = world_rank / 3; // Determine color based on row

    // Split the communicator based on the color and use the
    // original rank for ordering
    MPI_Comm row_comm;
    MPI_Comm_split(MPI_COMM_WORLD, color, world_rank, &row_comm);

    int row_rank, row_size;
    MPI_Comm_rank(row_comm, &row_rank);
    MPI_Comm_size(row_comm, &row_size);


    MPI_Comm copy_comm;
    MPI_Comm_dup(MPI_COMM_WORLD, &copy_comm);
    MPI_Group world_group;
    MPI_Comm_group(copy_comm, &world_group);

    int n = 4;
    const int ranks[4] = {0, 3, 6, 9};

     MPI_Group new_group;
     MPI_Group_incl(world_group, 4, ranks, &new_group);

    // Create new communicator based on group

     MPI_Comm new_comm;
     MPI_Comm_create(copy_comm, new_group, &new_comm);



    // Get the rank and size in the new communicator


    int new_rank = -1, new_size = -1;
    if (MPI_COMM_NULL != new_comm){
        MPI_Comm_rank(new_comm, &new_rank);
        MPI_Comm_size(new_comm, &new_size);
    }






    printf("WORLD RANK/SIZE: %d/%d \t ROW RANK/SIZE: %d/%d\t NEW RANK/SIZE: %d/%d\t color %d \n",
           world_rank, world_size, row_rank, row_size,  new_rank, new_size, color);

    MPI_Comm_free(&row_comm);
}
