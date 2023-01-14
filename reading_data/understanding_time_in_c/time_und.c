#include<stdio.h>

// mpicc -g -Wall -o mpi_peitro.out time_und.c
#include <sys/time.h>
#include <stdio.h>
 /*
 The struct timeval structure represents a calendar time. It has two members:
tv_sec : It is the number of seconds since the epoch.
tv_usec :It is additional microseconds after number of seconds calculation since the epoch. .
 */
int main() {
 
  struct timeval start, end;
  gettimeofday(&start, NULL);/*The 2nd argument points to the timezone structure. 
  It should normally be set to NULL because struct timezone is obsolete*/
 
  for (int i = 0; i <1e5 ; i++) {
  }
 
  gettimeofday(&end, NULL);//On success, the gettimeofday() return 0, for failure the function returns -1.
  printf("Time taken to count to 10^5 is : %ld micro seconds\n",
    ((end.tv_sec * 1000000 + end.tv_usec) -
    (start.tv_sec * 1000000 + start.tv_usec)));

  return 0;
}
