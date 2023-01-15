#include<stdio.h>

// mpicc -g -Wall -o mpi_peitro.out time_und.c
#include <sys/time.h>
#include <stdio.h>
 /*
 The struct timeval structure represents a calendar time. It has two members:
tv_sec : It is the number of seconds since the epoch.
tv_usec :It is additional microseconds after number of seconds calculation since the epoch. .
 */

/* function to measure time */
long int time_diff(struct timeval *start, struct timeval *end)
{
    return 1e+6 *(end->tv_sec - start->tv_sec) +  (end->tv_usec - start->tv_usec);
}
void convert_time_hour_sec(double seconds,long int *h,long int *m,long int *s)
{ 	
  *h = ((long int) seconds/3600); 
	*m = ((long int) seconds -(3600 * ( *h )))/60;
  *s = ((long int) seconds -(3600 * ( *h ))-(( *m ) * 60));
}
int main() {
 
  // struct timeval start, end;
  // gettimeofday(&start, NULL);/*The 2nd argument points to the timezone structure. 
  // It should normally be set to NULL because struct timezone is obsolete*/
 
  // for (int i = 0; i <1e5 ; i++) {
  // }
 
  // gettimeofday(&end, NULL);//On success, the gettimeofday() return 0, for failure the function returns -1.
  // printf("Time taken to count to 10^5 is : %ld micro seconds\n",
  //   (time_diff(&start,&end)));
  double t_nc_reading_time_Totalsum = 6000;
  long int t_seconds = 0;
  long int t_minutes = 0;
  long int t_hours = 0;
  convert_time_hour_sec(t_nc_reading_time_Totalsum,&t_hours,&t_minutes,&t_seconds);
  printf("The time taken to do Nc read is %ld hours,%ld minutes,%ld seconds \n", t_hours,t_minutes,t_seconds);
  return 0;
}
