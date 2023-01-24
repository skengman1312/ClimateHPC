#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mpi.h>
#include <omp.h>
#define LEVELS 5
#define GRID_POINTS 10
#define SUBARRAY_SIZE 2
/*
mpicc -std=c99 -g -Wall -fopenmp -I /apps/netCDF4.7.0--gcc-9.1.0/include -L /apps/netCDF4.7.0--gcc-9.1.0/lib -lnetcdf -o playing.out playing.c -lm 
mpirun -np 2 -prepend-rank ./playing.out

*/
int malloc2D(float ***array, int n, int m) {
    int i;
    /* allocate the n*m contiguous items */
    float *p = malloc(n*m*sizeof(float));
    if (!p) return -1;

    /* allocate the row pointers into the memory */
    (*array) = malloc(n*sizeof(float*));
    if (!(*array)) {
       free(p);
       return -1;
    }

    /* set up the pointers into the contiguous memory */
    for (i=0; i<n; i++)
       (*array)[i] = &(p[i*m]);

    return 0;
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int i,j;

    float **full_matrix;
    malloc2D(&full_matrix, LEVELS, GRID_POINTS); 
    int levels_per_proc = ceil((double)LEVELS / size);/*if 2.3 is passed to ceil(), it will return 3*/
    // if i have 2 process so i have 3 levels per process so i have 3*grixpoints on each process
    float **levels;
    malloc2D(&levels, levels_per_proc, GRID_POINTS);
    if (rank ==0){
        printf("\n");
        for (i = 0; i < LEVELS; i++){
            for (j = 0; j < GRID_POINTS; j++){
                full_matrix[i][j] = i*GRID_POINTS+j;
                printf("%f,",full_matrix[i][j]);
            }
        printf("\n");
        }
        printf("\n");
    }
    // int sizes[2] = {LEVELS, GRID_POINTS};
    // int subsizes[2] = {levels_per_proc, GRID_POINTS};
    // int starts[2] = {0, 0};
    int sendcounts[size];// for scatterv
    int displs[size];// for scatterv
    // MPI_Datatype subarray;
    // MPI_Type_create_subarray(SUBARRAY_SIZE,sizes,subsizes,starts,MPI_ORDER_C,MPI_FLOAT,&subarray);
    // // MPI_Type_create_resized(subarray, 0, GRID_POINTS*sizeof(float), &resizedtype);
    // MPI_Type_commit(&subarray);
       // intializing the arrays
    for (i = 0; i < size;i++){
        sendcounts[i]= levels_per_proc * GRID_POINTS;
        displs[i] = i * levels_per_proc * GRID_POINTS;
    }
    sendcounts[size-1] = (LEVELS-(levels_per_proc*(size-1)))*GRID_POINTS;
    

    MPI_Scatterv(&(full_matrix[0][0]),sendcounts,displs,MPI_FLOAT,&(levels[0][0]),GRID_POINTS*levels_per_proc,MPI_FLOAT,0,MPI_COMM_WORLD);
    if (rank ==0){
        for (i = 0; i < levels_per_proc; i++){
            for (j = 0; j < GRID_POINTS; j++){
                printf("%f,",levels[i][j]);
            }
            printf("\n");
        }
        printf("\n");
    }
    printf("\n");
    if (rank ==1){
        for (i = 0; i < levels_per_proc; i++){
            for (j = 0; j < GRID_POINTS; j++){
                printf("%f,",levels[i][j]);
            }
            printf("\n");
        }
        printf("\n");
    }

    MPI_Finalize();
    return 0;
}