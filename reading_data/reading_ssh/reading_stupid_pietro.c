//
// Created by Pietro on 17/01/2023.
//

#include <stdio.h>
#include <string.h>
#include <netcdf.h>
#include<mpi.h>
#include <string.h> // needed for memcpy

void get2DArray(int twoDArray[5]){
    // code to fill 2D array
}
void stupid_pietro(int twoDarray[5], int factor);


 int main(){
     // Declare the 3-dimensional array
     int threeDArray[4][5];
     // Declare the 2-dimensional array returned by the function
     int twoDArray[5];
     stupid_pietro(twoDArray, 1);
     // use memcpy to store 2DArray into the first slice of 3DArray
     memcpy(threeDArray[0], twoDArray, sizeof(twoDArray));
     for (int i = 0; i < 5; ++i) {
         printf("element of both arrays: %d, %d\n", threeDArray[0][i], twoDArray[i]);

     }

     stupid_pietro(twoDArray, 2);
     // use memcpy to store 2DArray into the first slice of 3DArray
     memcpy(threeDArray[1], twoDArray, sizeof(twoDArray));
     for (int i = 0; i < 5; ++i) {
         printf("element of both arrays: %d, %d\n", threeDArray[1][i], twoDArray[i]);

     }
     for (int i = 0; i < 2; ++i) {
         for (int j = 0; j < 5; ++j) {
             printf("element of 3d array: %d\n", threeDArray[i][j]);
         }
     }
     return 0;

 }

 void stupid_pietro(int twoDarray[5], int factor){
     for (int i = 0; i < 5; ++i) {
         twoDarray[i] = i*factor;

     }
};
