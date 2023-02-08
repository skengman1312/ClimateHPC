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
    for (int i = 0; i < 16; i) {
        printf("%d\t",i);
        l[j] = i;
        printf("%d\n", l[j]);
        j++;
        i += 4;


    }

    return 0;
};
