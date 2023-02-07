//
// Created by Pietro on 07/02/2023.
//
#include <stdio.h>
#include <string.h>
#include <netcdf.h>
#include <sys/time.h>
#include <stdlib.h>

int main(){
    int l[4] = {0};
    int j = 0;
    int k = 0;
    for (int i = 0; i < 12; i++) {
        printf("%d\n",i);
        k = i;
        if (i % 3 == 0) {
            l[j] = k;
            j++;
            printf("%d\t", l[j]);
        }
    }

    return 0;
};
