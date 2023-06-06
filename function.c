/**
 * @file: function.c
 * @author: Chanho Eom
 * @date: 8-Apr-2023
 * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <mpi.h>

#include "poisson1d.h"
#include "decomp1d.h"

void init_full_grid(double g[][maxn]){

int i,j;
const double junkval = -5;

for(i=0; i < maxn; i++){
for(j=0; j<maxn; j++){
g[i][j] = junkval;}}

}

void init_full_grids(double a[][maxn], double b[][maxn] ,double f[][maxn]){
 
int i,j;
const double junkval = -5;
  
for(i=0; i < maxn; i++){
for(j=0; j<maxn; j++){
a[i][j] = junkval;
b[i][j] = junkval;
f[i][j] = junkval;}}

}


/* prints to stdout in GRID view */
void print_full_grid(double x[][maxn], int nx){

int i,j;
 
for(j = nx + 1; j>=0; j--){
for(i = 0; i < nx +2; i++){
if(x[i][j] < 10000.0){
printf("|%2.6lf| ",x[i][j]);
} else {
printf("%9.2lf ",x[i][j]);}}
printf("\n");}

}

void init_basic_1d(double a[][maxn], double b[][maxn], double f[][maxn], int nx, int ny, int s, int e){

//printf("================== Initialising Poisson's boundary conditions starts ==================\n");

double x, y;

for (int i = s -1; i <= e +1; i++){
for (int j = 0; j <= nx +1; j++){
a[i][j] = 0.0;
b[i][j] = 0.0;
f[i][j] = 0.0;}}

for (int i = s; i <= e; i++){
x = 1.0/((double)nx + 1.0) * i;
a[i][0] = 0.0;
b[i][0] = 0.0;
a[i][nx + 1] = 1.0/((x + 1.0)*(x + 1.0) + 1.0);
b[i][nx + 1] = 1.0/((x + 1.0)*(x + 1.0) + 1.0);
}

if (s == 1){
for (int j = 1; j < nx +2; j++){
y = 1.0/((double) ny + 1.0) * j;
a[0][j] = y / (y*y + 1.0);
b[0][j] = y / (y*y + 1.0);
}}

if (e == nx){
for (int j = 1; j < nx +2; j++){
y = 1.0/((double)ny + 1.0) * j;
a[nx + 1][j] = y / (y*y + 4.0);
b[nx + 1][j] = y / (y*y + 4.0);
}} 

}

void init_basic_2d(double a[][maxn], double b[][maxn], double f[][maxn], int nx, int ny, int s[2], int e[2]){

double x, y;

for (int i = s[0] -1; i <= e[0] +1; i++){
for (int j = s[1] -1; j <= e[1] +1; j++){
a[i][j] = 0.0;
b[i][j] = 0.0;
f[i][j] = 0.0;}}

if (s[0] == 1){
for (int j = s[1] -1; j <= e[1] +1; j++){
y = 1.0 / ((double)ny + 1.0) * j;
a[0][j] = y / (y*y + 1.0);
b[0][j] = y / (y*y + 1.0);}
}

if (s[1] == 1){
for (int i = s[0] -1; i<= e[0] +1; i++){
a[i][0] = 0.0;
b[i][0] = 0.0;}
}

if (e[0] == nx){
for (int j = s[1] -1; j<= e[1] +1; j++){
y = 1.0 / ((double)ny + 1.0) * j;
a[nx +1][j] = y / (4 + y*y);
b[nx +1][j] = y / (4 + y*y);}
}

if (e[1] == ny){
for (int i = s[0] -1; i <= e[0] +1; i++){ 
x = 1.0 / ((double)nx + 1.0) * i;
a[i][ny +1] = 1 / ((x +1)*(x +1) + 1);
b[i][ny +1] = 1 / ((x +1)*(x +1) + 1);}
}

}


void onedinit_basic(double a[][maxn], double b[][maxn], double f[][maxn], int nx, int ny, int s, int e){

//printf("================== Initialising Given boundary conditions starts ==================\n");

int i,j;
double left, bottom, right, top;

left   = -1.0;
bottom = 1.0;
right  = 2.0;
top    = 3.0;

/* set everything to 0 first */
for(i=s-1; i<=e+1; i++){
for(j=0; j <= nx+1; j++){
      a[i][j] = 0.0;
      b[i][j] = 0.0;
      f[i][j] = 0.0;}}

  /* deal with boundaries */
for(i=s; i<=e; i++){
    a[i][0] = bottom;
    b[i][0] = bottom;
    a[i][nx+1] = top;
    b[i][nx+1] = top;}

  /* this is true for proc 0 */
if( s == 1 ){
for(j=1; j<nx+1; j++){
      a[0][j] = left;
      b[0][j] = left;}}

  /* this is true for proc size-1 */
if( e == nx ){
for(j=1; j<nx+1; j++){
      a[nx+1][j] = right;
      b[nx+1][j] = right;}}

}

void print_in_order(double x[][maxn], MPI_Comm comm, int nx){

//printf("================== Printing in order starts ==================\n");

int myid, size;
int i;

  
MPI_Comm_rank(comm, &myid);
MPI_Comm_size(comm, &size);  

MPI_Barrier(comm); 
printf("Attempting to print in order\n");  
sleep(1);
  
MPI_Barrier(comm);
  
for(i=0; i<size; i++){
if( i == myid ){
printf("\n=================== proc %d ==================\n",myid);
print_full_grid(x, nx);}
fflush(stdout);
usleep(500);
MPI_Barrier(MPI_COMM_WORLD);}

}

void print_grid_to_file(char *fname, double x[][maxn], int nx, int ny){

//printf("================== Printing the grid to the file starts ==================\n");

FILE *fp;
int i,j;

fp = fopen(fname, "w");
if( !fp ){
fprintf(stderr, "Error: can't open file %s\n",fname);
exit(4);}

for(j=ny+1; j>=0; j--){
for(i=0; i<nx+2; i++){
fprintf(fp, "%lf ",x[i][j]);}
fprintf(fp, "\n");}
fclose(fp);

}

void solution_grid(double solution[][maxn], int nx, int ny){

double x, y;
for( int i = 0; i < nx + 2; i++){
for (int j = 0; j < ny + 2; j++){
x = 1.0 / ((double)nx + 1.0) *i;
y = 1.0 / ((double)ny + 1.0) * j;
solution[i][j] = y / ((1.0 + x) * (1.0 + x) + y * y);}}

}

void gather_grid(double a[][maxn], int nx, int nprocs, int myid, MPI_Comm comm){

int s[12], e[12], cols[12];	// Chuck is able to use  at most 12 processors

if (myid == 0){
if (nprocs > 12){
perror("The number of processors is more than 12");
exit(EXIT_FAILURE);}}

for (int i = 0; i < nprocs; i++){
decomp_1d(nx, i, nprocs, (s + i), (e + i));
cols[i] = e[i] - s[i] +1;}

if (myid != 0){
for (int j = 0; j < cols[myid]; j++){
MPI_Send(&a[s[myid] +j][0], (nx +2), MPI_DOUBLE, 0, myid * j + myid, comm);
}}

if (myid == 0){
for (int i = 1; i < nprocs; i++){
for (int j = 0; j < cols[i]; j++){
MPI_Recv(&a[s[i] +j][0], (nx +2), MPI_DOUBLE, i, i*j + i, comm, MPI_STATUS_IGNORE);
}}}

// For the right boundary
if (myid == nprocs -1){
MPI_Send(&a[nx +1], (nx +2), MPI_DOUBLE, 0, 99, comm);}
if (myid == 0){
MPI_Recv(&a[nx +1], (nx +2), MPI_DOUBLE, nprocs -1, 99, comm, MPI_STATUS_IGNORE);} 

}


void gather_grid_2d(double a[][maxn], int nx, int ny, int nprocs, int myid, int s[2], int e[2], MPI_Comm comm){

int s_proc[16][2], e_proc[16][2], cols, rows, cols_proc[16], rows_proc[16];

if (myid == 0){
if (nprocs > 16){
perror("The number of processors is more than 12");
exit(EXIT_FAILURE);}}

rows = e[1] - s[1] +1;

if (myid == 0){
for (int j = 0; j < 2; j++){
s_proc[0][j] = s[j];
e_proc[0][j] = e[j];}

for (int i = 1; i < nprocs; i++){
MPI_Recv(s_proc[i], 2, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
MPI_Recv(e_proc[i], 2, MPI_INT, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);}

for (int i = 0; i < nprocs; i++){
cols_proc[i] = e_proc[i][0] - s_proc[i][0] +1;
rows_proc[i] = e_proc[i][1] - s_proc[i][1] +1;
printf("(myid: %d) x_axis: (s[%d][0] = %d, e[%d][0] = %d), y_axis: (s[%d][1] = %d, e[%d][1] = %d), num_cols[%d] = %d, num_rows[%d] = %d \n", i, i, s_proc[i][0], i, e_proc[i][0], i, s_proc[i][1], i, e_proc[i][1], i, cols_proc[i], i, rows_proc[i]);}

}

if (myid != 0){
MPI_Send(s, 2, MPI_INT, 0, 0, MPI_COMM_WORLD);
MPI_Send(e, 2, MPI_INT, 0, 1, MPI_COMM_WORLD);}

if (myid != 0){

double sendbuf[maxn * maxn];
cols = e[0] - s[0] +1;

for (int i = s[0]; i <= e[0]; i++){
for (int j = s[1]; j <= e[1]; j++){
sendbuf[(i - s[0]) * rows + (j - s[1])] = a[i][j];}}

MPI_Send(sendbuf, cols * rows, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
}

if (myid == 0){

double recvbuf[maxn * maxn];

for (int i = 1; i < nprocs; i++){
MPI_Recv(recvbuf, cols_proc[i] * rows_proc[i], MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

for (int j = s_proc[i][0]; j <= e_proc[i][0]; j++){
for (int k = s_proc[i][1]; k <= e_proc[i][1]; k++){
a[j][k] = recvbuf[(j - s_proc[i][0]) * rows_proc[i] + (k - s_proc[i][1])];}}
}}



// Cheating in a sense to get the boundary condition from the analytic solution
double b[maxn][maxn];
solution_grid(b, nx, ny);

if (myid == 0){
for (int i = 0; i <= nx +1; i++){
a[i][0] = b[i][0];
a[0][i] = b[0][i];
a[nx +1][i] = b[nx +1][i];
a[i][nx +1] = b[i][nx +1];}}

}


void write_grid(double a[][maxn], int nx, int myid, int s, int e){

FILE *fp;
char fname[50] = "Grid_rank";
char rank[7];

sprintf(rank, "%d", myid);
strcat(fname, rank);
strcat(fname, ".txt");

fp = fopen(fname, "w");
if ( !fp ){
fprintf(stderr, "Error for opening the file %s\n", fname);
exit(4); }

if ( s == 1){
for (int j = nx +1; j >= 0; j--){
for (int i = 0; i <= e; i++){
fprintf(fp, "%lf", a[i][j]);}
fprintf(fp, "\n");}

fclose(fp);
} else if (e == nx){
for (int j = nx +1; j >= 0; j--){
for (int i = s; i <= nx +1; i++){
fprintf(fp, "%lf", a[i][j]);}
fprintf(fp, "\n");}

fclose(fp);
} else{
for (int j = nx +1; j >= 0; j--){
for (int i = s; i <= e; i++){
fprintf(fp, "%lf", a[i][j]);}
fprintf(fp, "\n");}

fclose(fp);}

}


