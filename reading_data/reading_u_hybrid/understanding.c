#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>


int main(int argc, char *argv[])
{
    /* MPI  inizialization */
    MPI_Init(&argc, &argv);

    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int x[5][3] = {{1,2,3}, {4,5,6},{7,8,9},{10,11,12},{13,14,15}};
    int size2 = 5 / size;
    int y[size2][3] ;
    int sendcnt=size2*3;
    int recvcnt=size2*3;
    MPI_Scatter(x, sendcnt, MPI_INT,
                y, recvcnt, MPI_INT, 0, MPI_COMM_WORLD);
    int i = 0;
    for (i = 0; i < 3;i++)
    {
        printf("%d for process %d\n", y[0][i],rank);
    }
    MPI_Finalize();
    return 0;
}