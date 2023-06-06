/**
 * @file: decomp1d.h
 * @author: Chanho Eom
 * @date: 8-Apr-2023
 * */

#ifndef DECOMP1D_H
#define DECOMP1D_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#include "poisson1d.h"

void decomp_1d(int nx, int myid, int nprocs, int *s, int *e);
void decomp_2d(int nx, int myid, int dims[2], int coords[2], int s[2], int e[2]);



#endif 
