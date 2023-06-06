/**
 * @file: decomp1d.c
 * @author: Chanho Eom
 * @date: 8-Apr-2023
 * */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#include "poisson1d.h"
#include "decomp1d.h"

void decomp_1d(int nx, int myid, int nprocs, int *s, int *e){

int size = nx / nprocs;
int lest = nx % nprocs;

if (myid < lest){
*s = myid * size + 1 + myid;
*e = *s + size;
} else{
*s = myid * size + 1 + lest;
*e = *s + (size - 1);
}

if (*e > nx || myid == nprocs - 1){
*e = nx;
}

}

void decomp_2d(int nx, int myid, int dims[2], int coords[2], int s[2], int e[2]){

int rank_1d = (myid - coords[1]) / dims[1];
decomp_1d(nx, rank_1d, dims[0], s, e);

rank_1d = (myid) % dims[1];

decomp_1d(nx, rank_1d, dims[1],  &s[1], &e[1]);

}

