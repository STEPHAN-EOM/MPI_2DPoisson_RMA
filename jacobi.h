/**
 * @file: jacobi.h
 * @author: Chanho Eom
 * @date: 23-Apr-2023
 * */

#ifndef JACOBI_H
#define JACOBI_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#include "poisson1d.h"



void sweep1d(double a[][maxn], double f[][maxn], int ny, int s, int e, double b[][maxn]);

void exchang1(double x[][maxn], int nx, int s, int e, MPI_Comm comm, int nbrbottom, int nbrtop);
void exchang2(double x[][maxn], int nx, int s, int e, MPI_Comm comm, int nbrleft, int nbrright);
void exchang3(double x[][maxn], int nx, int s, int e, MPI_Comm comm, int nbrleft, int nbrright);
void exchangi1(double x[][maxn], int nx, int s, int e, MPI_Comm comm, int nbrleft, int nbrright);

// For the 2D Ghost Exchange Routine with RMA by C. Eom
void exchang_2d_RMA(double x[][maxn], int nx, int s[2], int e[2], MPI_Comm comm, int nbrleft, int nbrright, int nbrdown, int nbrup, MPI_Win win);

// For the 2D Ghost Exchange Routine with RMA PSCW by C. Eom
void exchang_2d_RMA_pscw(double x[][maxn], int nx, int s[2], int e[2], MPI_Comm comm, int nbrleft, int nbrright, int nbrdown, int nbrup, MPI_Win win, MPI_Group origin[2], MPI_Group target[2]);


// For 2d Ghost Exchange Routine with Sendrecv by C. Eom
void exchang3_2d(double x[][maxn], int nx, int s[2], int e[2], MPI_Comm comm, int nbrleft, int nbrright, int nbrdown, int nbrup);
void sweep_2d(double a[][maxn], double f[][maxn], int nx, int s[2], int e[2], double b[][maxn]);


// For 2d Ghost Exchange Routine with Non-blocking send and recv by C. Eom
void nbxchange_and_sweep_2d(double u[][maxn], double f[][maxn], int nx, int ny, int s[2], int e[2], double unew[][maxn], MPI_Comm comm, int nbrleft, int nbrright, int nbrdown, int nbrup);


void nbxchange_and_sweep(double u[][maxn], double f[][maxn], int nx, int ny,
			 int s, int e, double unew[][maxn], MPI_Comm comm,
			 int nbrleft, int nbrright);

double griddiff(double a[][maxn], double b[][maxn], int nx, int s, int e);

// For 2d Decomposition by C. Eom
double griddiff_2d(double a[][maxn], double b[][maxn], int nx, int s[2], int e[2]);

#endif
