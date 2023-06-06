/**
 * @file: function.h
 * @author: Chanho Eom
 * @date: 8-Apr-2023
 * */

#ifndef FUNCTION_H
#define FUNCTION_H


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <mpi.h>

#include "poisson1d.h"
#include "decomp1d.h"

void init_full_grid(double g[][maxn]);
void init_full_grids(double a[][maxn], double b[][maxn], double f[][maxn]); 

void init_basic_1d(double a[][maxn], double b[][maxn], double f[][maxn], int nx, int ny, int s, int e);
void init_basic_2d(double a[][maxn], double b[][maxn], double f[][maxn], int nx, int ny, int s[2], int e[2]);
void onedinit_basic(double a[][maxn], double b[][maxn], double f[][maxn],int nx, int ny, int s, int e);

void print_full_grid(double x[][maxn], int nx);
void print_in_order(double x[][maxn], MPI_Comm comm, int nx);
void print_grid_to_file(char *fname, double x[][maxn], int nx, int ny);

void solution_grid(double solution[][maxn], int nx, int ny);
void gather_grid(double a[][maxn], int nx, int nprocs, int myid, MPI_Comm comm);
void gather_grid_2d(double a[][maxn], int nx, int ny, int nprocs, int myid, int s[2], int e[2], MPI_Comm comm);
void write_grid(double a[][maxn], int nx, int s, int e, int myid);

#endif

