#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int global_matrix[4][6] = { {1, 2, 3, 4, 5, 6},
                                {7, 8, 9, 10, 11, 12},
                                {13, 14, 15, 16, 17, 18},
                                {19, 20, 21, 22, 23, 24} };

    int local_matrix[2][6];
    int counts[size], displs[size];

    // Set the counts and displacements for each process
    for (int i = 0; i < size; i++) {
        counts[i] = 2;
        displs[i] = i * 2;
    }

    // Create the subarray datatype
    MPI_Datatype subarray;
    int sizes[2] = {4, 6};
    int subsizes[2] = {2, 6};
    int starts[2] = {0, 0};
    MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_C, MPI_INT, &subarray);
    MPI_Type_commit(&subarray);

    // Scatter the global matrix to the local matrices
    MPI_Scatterv(global_matrix, counts, displs, subarray, local_matrix, 2*6, MPI_INT, 0, MPI_COMM_WORLD);

    // Print the local matrix on each process
    printf("Process %d has the following local matrix: \n", rank);
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 6; j++) {
            printf("%d ", local_matrix[i][j]);
        }
        printf("\n");
    }

    // Clean up
    MPI_Type_free(&subarray);
    MPI_Finalize();
    return 0;
}