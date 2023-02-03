#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define BIGSIZE 5
#define GRID_POINTS 3
int malloc2D(float ***array, int n, int m);

int main(int argc, char *argv[])
{
    /* MPI  inizialization */
    MPI_Init(&argc, &argv);

    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // int i, j;
    // int x[5][3] = {{1,2,3}, {4,5,6},{7,8,9},{10,11,12},{13,14,15}};
    // float **x = (float**) malloc(sizeof(float *) * BIGSIZE);
    // for (int i = 0; i < BIGSIZE; i++){
    //     x[i] = (float *)calloc(GRID_POINTS ,sizeof(float));
    // }
    float **x;
    malloc2D(&x,BIGSIZE,GRID_POINTS);
    /*INtialize */
    if (rank ==0){
        for (int i = 0; i < BIGSIZE; i++){
            for (int j = 0; j < GRID_POINTS; j++){
                x[i][j] = i*GRID_POINTS+j;
                printf("%f,",x[i][j]);
            }
                    printf("\n");
        }
        printf("\n");
    }
    /* you have to ask for 5 processes*/
    float y[GRID_POINTS] ;
    int sendcnt=GRID_POINTS;
    int recvcnt=GRID_POINTS;
    MPI_Scatter(&x[0][0], sendcnt, MPI_INT,
                y, recvcnt, MPI_INT, 0, MPI_COMM_WORLD);
    for (int i = 0; i < 3; i++)
    {
        printf("%f for process %d\n", y[i],rank);
    }
    MPI_Finalize();
    return 0;
}
int malloc2D(float ***array, int n, int m) {
    int i,local_ok;
    /* allocate the n*m contiguous items */
    float *p = calloc(n*m,sizeof(float));

    if (!p) return -1;

    /* allocate the row pointers into the memory */
    (*array) = malloc(n*sizeof(float*));
    if (!(*array)) {
       free(p);
       return -1;
    }
    // if ( ( p == NULL ) || ( array == NULL ) ) {
    //         local_ok = 0;
            
    // }
    // Check_for_error(local_ok, "Main for loop", "Can't allocate memory vector", MPI_COMM_WORLD);

    /* set up the pointers into the contiguous memory */
    for (i=0; i<n; i++)
       (*array)[i] = &(p[i*m]);

    return 0;
}