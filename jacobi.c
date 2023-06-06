/**
 * @file: jacobi.c
 * @author: Chanho Eom
 * @date: 23-Apr-2023
 * */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>

#include "poisson1d.h"
#include "jacobi.h"

/* sequentialized if there is no buffering */
void exchang1(double x[][maxn], int nx, int s, int e, MPI_Comm comm, int nbrleft, int nbrright)
{
  int rank;
  int ny;

  MPI_Comm_rank(comm, &rank);

  ny = nx;

  MPI_Ssend(&x[e][1], ny, MPI_DOUBLE, nbrright, 0, comm);
  /* printf("(myid: %d) sent \"col\" %d with %d entries to nbr: %d\n",rank, e, ny, nbrright); */
  
  MPI_Recv(&x[s-1][1], ny, MPI_DOUBLE, nbrleft, 0, comm, MPI_STATUS_IGNORE);
  /* printf("(myid: %d) recvd into \"col\" %d from %d entries from nbr: %d\n",rank, s-1, ny, nbrleft); */

  MPI_Ssend(&x[s][1], ny, MPI_DOUBLE, nbrleft, 1, comm);
  /* printf("(myid: %d) sent \"col\" %d with %d entries to nbr: %d\n",rank, s, ny, nbrleft); */
  MPI_Recv(&x[e+1][1], ny, MPI_DOUBLE, nbrright, 1, comm, MPI_STATUS_IGNORE);
  /* printf("(myid: %d) recvd into \"col\" %d from %d entries from nbr: %d\n",rank, e+1, ny, nbrright); */

}

void sweep1d(double a[][maxn], double f[][maxn], int nx, int s, int e, double b[][maxn])
{
  double h;
  int i,j;

  h = 1.0/((double)(nx+1));

  for(i=s; i<=e; i++){
    for(j=1; j<nx+1; j++){
      b[i][j] = 0.25 * ( a[i-1][j] + a[i+1][j] + a[i][j+1] + a[i][j-1]  - h*h*f[i][j] );
    }
  }
}

/* ordered sends / receives */
void exchang2(double x[][maxn], int nx, int s, int e, MPI_Comm comm, int nbrleft, int nbrright)
{
  int coord;
  int rank;

  MPI_Comm_rank(comm, &rank);

  coord = rank;

  if(coord%2 == 0){

    MPI_Ssend(&x[e][1], nx, MPI_DOUBLE, nbrright, 0, comm);

    MPI_Recv(&x[s-1][1], nx, MPI_DOUBLE, nbrleft, 0, comm, MPI_STATUS_IGNORE);

    MPI_Ssend(&x[s][1], nx, MPI_DOUBLE, nbrleft, 1, comm);

    MPI_Recv(&x[e+1][1], nx, MPI_DOUBLE, nbrright, 1, comm, MPI_STATUS_IGNORE);
  } else {
    MPI_Recv(&x[s-1][1], nx, MPI_DOUBLE, nbrleft, 0, comm, MPI_STATUS_IGNORE);

    MPI_Ssend(&x[e][1], nx, MPI_DOUBLE, nbrright, 0, comm);

    MPI_Recv(&x[e+1][1], nx, MPI_DOUBLE, nbrright, 1, comm, MPI_STATUS_IGNORE);

    MPI_Ssend(&x[s][1], nx, MPI_DOUBLE, nbrleft, 1, comm);

  }

}

/* sendrecv */
void exchang3(double x[][maxn], int nx, int s, int e, MPI_Comm comm, int nbrleft, int nbrright)
{

  MPI_Sendrecv(&x[e][1], nx, MPI_DOUBLE, nbrright, 0, &x[s-1][1], nx, MPI_DOUBLE, nbrleft,
	       0, comm, MPI_STATUS_IGNORE);

  MPI_Sendrecv(&x[s][1], nx, MPI_DOUBLE, nbrleft, 1, &x[e+1][1], nx, MPI_DOUBLE, nbrright,
	       1, comm, MPI_STATUS_IGNORE);

}

void exchangi1(double x[][maxn], int nx, int s, int e, MPI_Comm comm, int nbrleft, int nbrright)
{
  MPI_Request reqs[4];

  MPI_Irecv(&x[s-1][1], nx, MPI_DOUBLE, nbrleft, 0, comm, &reqs[0]);
  MPI_Irecv(&x[e+1][1], nx, MPI_DOUBLE, nbrright, 0, comm, &reqs[1]);
  MPI_Isend(&x[e][1],   nx, MPI_DOUBLE, nbrright, 0, comm, &reqs[2]);
  MPI_Isend(&x[s][1],   nx, MPI_DOUBLE, nbrleft, 0, comm, &reqs[3]);
  /* not doing anything useful here */

  MPI_Waitall(4, reqs, MPI_STATUSES_IGNORE);
}

// Exercise 2. Question.1
void exchang_2d_RMA(double x[][maxn], int nx, int s[2], int e[2], MPI_Comm comm, int nbrleft, int nbrright, int nbrdown, int nbrup, MPI_Win win){

// Define a vector for the row exchanges
MPI_Datatype newtype;
MPI_Type_vector(e[0] - s[0] +1, 1, maxn, MPI_DOUBLE, &newtype);
MPI_Type_commit(&newtype);

// Integer type that can hold an arbitrary memory address; Displacement from start of window to target buffer
MPI_Aint offset;
MPI_Win_fence(MPI_MODE_NOPRECEDE, win);

// left direction
offset = s[1];
MPI_Put(&x[e[0]][s[1]], e[1] - s[1] +1, MPI_DOUBLE, nbrright, offset, e[1] - s[1] +1, MPI_DOUBLE, win);

// down direction
offset = maxn + e[1];
//offset = s[0];
MPI_Put(&x[s[0]][e[1]], 1, newtype, nbrup, offset, 1, newtype, win);

// right direction
offset = maxn + s[1];
//offset = s[0];
//MPI_Get(&x[s[0]][s[1]], e[1] - s[1] +1, MPI_DOUBLE, nbrleft, offset, e[1] - s[1] +1, MPI_DOUBLE, win);
MPI_Get(&x[e[0] +1][s[1]], e[1] - s[1] +1, MPI_DOUBLE, nbrright, offset, e[1] - s[1] +1, MPI_DOUBLE, win);

// up direction
offset = maxn + e[1] +1;
//offset = e[0];
//MPI_Get(&x[s[0]][s[1]], 1, newtype, nbrdown, offset, 1, newtype, win);
MPI_Get(&x[s[0]][e[1] +1], 1, newtype, nbrup, offset, 1, newtype, win);

MPI_Win_fence(MPI_MODE_NOPRECEDE, win);

MPI_Type_free(&newtype);
}


// Exercise 2. Question.2 
void exchang_2d_RMA_pscw(double x[][maxn], int nx, int s[2], int e[2], MPI_Comm comm, int nbrleft, int nbrright, int nbrdown, int nbrup, MPI_Win win, MPI_Group origin[2], MPI_Group target[2]){

MPI_Datatype newtype;
MPI_Type_vector(e[0] - s[0] +1, 1, maxn, MPI_DOUBLE, &newtype);
MPI_Type_commit(&newtype);

MPI_Aint offset;

if (nbrleft != MPI_PROC_NULL){
MPI_Win_post(target[0], 0, win);
MPI_Win_wait(win);}

if (nbrdown != MPI_PROC_NULL){
MPI_Win_post(target[1], 0, win);
MPI_Win_wait(win);}

if (nbrright != MPI_PROC_NULL){
MPI_Win_start(origin[0], 0, win);

// left direction
offset = s[1];
MPI_Put(&x[e[0]][s[1]], e[1] - s[1] +1, MPI_DOUBLE, nbrright, offset, e[1] - s[1] +1, MPI_DOUBLE, win);

// right direction
offset = maxn + s[1];
//offset = s[0];
//MPI_Get(&x[s[0]][s[1]], e[1] - s[1] +1, MPI_DOUBLE, nbrleft, offset, e[1] - s[1] +1, MPI_DOUBLE, win);
MPI_Get(&x[e[0] +1][s[1]], e[1] - s[1] +1, MPI_DOUBLE, nbrleft, offset, e[1] - s[1] +1, MPI_DOUBLE, win); 
MPI_Win_complete(win);
}

if (nbrup != MPI_PROC_NULL){
MPI_Win_start(origin[1], 0, win);

// down direction
offset = maxn + e[1];
//offset = s[0];
MPI_Put(&x[s[0]][e[1]], 1, newtype, nbrup, offset, 1, newtype, win);

// up direction
offset = maxn + e[1] +1;
// offset = e[0];
// MPI_Get(&x[s[0]][s[1]]
MPI_Get(&x[s[0]][e[1] +1], 1, newtype, nbrdown, offset, 1, newtype, win);
MPI_Win_complete(win);
}

}


void exchang3_2d(double x[][maxn], int nx, int s[2], int e[2], MPI_Comm comm, int nbrleft, int nbrright, int nbrdown, int nbrup){

MPI_Sendrecv(&x[e[0]][s[1]], e[1] - s[1] +1, MPI_DOUBLE, nbrright, 0, &x[s[0] -1][s[1]], e[1] - s[1] +1, MPI_DOUBLE, nbrleft, 0, comm, MPI_STATUS_IGNORE);
MPI_Sendrecv(&x[s[0]][s[1]], e[1] - s[1] +1, MPI_DOUBLE, nbrleft, 1, &x[e[0] +1][s[1]], e[1] - s[1] +1, MPI_DOUBLE, nbrright, 1, comm, MPI_STATUS_IGNORE);

MPI_Datatype newtype;
MPI_Type_vector(e[0] - s[0] +1, 1, maxn, MPI_DOUBLE, &newtype);
MPI_Type_commit(&newtype);

MPI_Sendrecv(&x[s[0]][s[1]], 1, newtype, nbrdown, 2, &x[s[0]][s[1] -1], 1, newtype, nbrup, 2, comm, MPI_STATUS_IGNORE);
MPI_Sendrecv(&x[s[0]][e[1]], 1, newtype, nbrup, 3, &x[s[0]][e[1] +1], 1, newtype, nbrdown, 3, comm, MPI_STATUS_IGNORE);

}

void sweep_2d(double a[][maxn], double f[][maxn], int nx, int s[2], int e[2], double b[][maxn]){
 
double h;
  
h = 1.0/((double)(nx+1));
  
for(int i = s[0]; i <= e[0]; i++){
for(int j = s[1]; j <= e[1]; j++){
b[i][j] = 0.25 * ( a[i-1][j] + a[i+1][j] + a[i][j+1] + a[i][j-1]  - h*h*f[i][j] );}}

}

void nbxchange_and_sweep_2d(double u[][maxn], double f[][maxn], int nx, int ny, int s[2], int e[2], double unew[][maxn], MPI_Comm comm, int nbrleft, int nbrright, int nbrdown, int nbrup){

MPI_Request req[8];
MPI_Status status;
int idx;
double h;

int myid;
MPI_Comm_rank(comm, &myid);

h = 1.0/( (double)(nx+1));

MPI_Irecv(&u[s[0] -1][s[1]], ny, MPI_DOUBLE, nbrleft, 1, comm, &req[0] );
MPI_Irecv(&u[e[0] +1][s[1]], ny, MPI_DOUBLE, nbrright, 2, comm, &req[1] );
MPI_Isend(&u[e[0]][s[1]], ny, MPI_DOUBLE, nbrright, 1, comm, &req[2]);
MPI_Isend(&u[s[0]][s[1]], ny, MPI_DOUBLE, nbrleft, 2, comm, &req[3]);


MPI_Datatype newtype;
MPI_Type_vector(e[0] - s[0] +1, 1, maxn, MPI_DOUBLE, &newtype);
MPI_Type_commit(&newtype);


MPI_Irecv(&u[s[0]][s[1] -1], 1, newtype, nbrdown, 3, comm, &req[4] );
MPI_Irecv(&u[s[0]][e[1] +1], 1, newtype, nbrup, 4, comm, &req[5] );
MPI_Isend(&u[s[0]][e[1]], 1, newtype, nbrup, 3, comm, &req[6]);
MPI_Isend(&u[s[0]][s[1]], 1, newtype, nbrdown, 4, comm, &req[7]);



/* perform purely local updates (that don't need ghosts) */
/* 2 cols or less means all are on processor boundary */
if((e[0] - s[0] +1 >= 3) && (e[1] - s[1] +1 >= 3)){
for(int i = s[0] +1; i <= e[0]; i++){
for(int j = s[1] +1; j < e[1]; j++){
unew[i][j] = 0.25 * ( u[i-1][j] + u[i+1][j] + u[i][j+1] + u[i][j-1]  - h*h*f[i][j] );}}
}


// For the left boundary that does not need ghost
if ((s[0] == 1) && (e[0] - s[0] +1 >= 2)){
for (int j = s[1] + 1; j < e[1]; j++){
unew[1][j] = 0.25 * ( u[0][j] + u[2][j] + u[1][j -1] + u[1][j +1] - h*h*f[1][j] );}}

// For the right boundary that does not need ghost
if ((e[0] == nx) && (e[0] - s[0] +1 >= 2)){
for (int j = s[1] + 1; j < e[1]; j++){
unew[nx][j] = 0.25 * ( u[nx -1][j] + u[nx +1][j] + u[nx][j -1] + u[nx][j +1] - h*h*f[nx][j] );}}

// For the bottom boundary that does not need ghost
if ((s[1] == 1) && (e[1] - s[1] +1 >= 2)){
for (int i = s[0] +1; i < e[1]; i++){
unew[i][1] = 0.25 * ( u[i -1][1] + u[i +1][1] + u[i][0] + u[i][2] - h*h*f[i][1] );}}

// For the top boundary that does not need ghost
if ((e[1] == ny) && (e[1] - s[1] +1 >= 2)){
for (int i = s[0] +1; i < e[0]; i++){
unew[i][ny] = 0.25 * ( u[i -1][ny] + u[i +1][ny] + u[i][ny -1] + u[i][ny +1] - h*h*f[i][ny] );}}


for(int k=0; k < 8; k++){
MPI_Waitany(8, req, &idx, &status);


// idx 0, 1 are recvs for left, right
switch(idx){
case 0:

if (nbrleft != MPI_PROC_NULL && (e[0] -s[0] +1 > 1)){
for (int j = s[1] +1; j < e[1]; j++){
unew[s[0]][j] = 0.25 * ( u[s[0] -1][j] + u[s[0] +1][j] + u[s[0]][j -1] + u[s[0]][j +1] - h*h*f[s[0]][j] );
}}
break;

case 1:

if (nbrright != MPI_PROC_NULL && (e[0] - s[0] +1 > 1)){
for(int j = s[1] +1; j < e[1]; j++){       
unew[e[0]][j] = 0.25 * ( u[e[0] -1][j] + u[e[0] +1][j] + u[e[0]][j -1] + u[e[0]][j +1] - h*h*f[e[0]][j] ); 
}}
break;

// idx 4, 5 are recvs for down, up
case 4:

if (nbrdown != MPI_PROC_NULL && (e[1] - s[1] +1 > 1)){
for (int i = s[0] +1; i < e[0]; i++){
unew[i][s[1]] = 0.25 * ( u[i -1][s[1]] + u[i +1][s[1]] + u[i][s[1] -1] + u[i][s[1] +1] - h*h*f[i][s[1]] );
}}
break;

case 5:

if (nbrup != MPI_PROC_NULL && (e[1] - s[1] +1 > 1)){
for (int i = s[0] +1; i < e[0]; i++){
unew[i][e[1]] = 0.25 * ( u[i -1][e[1]] + u[i +1][e[1]] + u[i][e[1] -1] + u[i][e[1] +1] - h*h*f[i][e[1]] ); 
}}
break;
      
default:      
break;
}}


// Update the corners. 
unew[s[0]][s[1]] = 0.25 * ( u[s[0] -1][s[1]] + u[s[0] +1][s[1]] + u[s[0]][s[1] -1] + u[s[0]][s[1] +1] - h*h*f[s[0]][s[1]] );	// (Start, Start)
unew[s[0]][e[1]] = 0.25 * ( u[s[0] -1][e[1]] + u[s[0] +1][e[1]] + u[s[0]][e[1] -1] + u[s[0]][e[1] +1] - h*h*f[s[0]][e[1]] );	// (Start, End)
unew[e[0]][s[1]] = 0.25 * ( u[e[0] -1][s[1]] + u[e[0] +1][s[1]] + u[e[0]][s[1] -1] + u[e[0]][s[1] +1] - h*h*f[e[0]][s[1]] );	// (End, Start)
unew[e[0]][e[1]] = 0.25 * ( u[e[0] -1][e[1]] + u[e[0] +1][e[1]] + u[e[0]][e[1] -1] + u[e[0]][e[1] +1] - h*h*f[e[0]][e[1]] );	// (End, End)

if (e[1] - s[1] +1 == 1) {
for (int i = s[0] +1; i < e[0]; i++) {
unew[i][s[1]] = 0.25 * ( u[i -1][s[1]] + u[i +1][s[1]] + u[i][s[1] -1] + u[i][s[1] +1] - h*h*f[i][s[1]] );
}}

if (e[0] -s[0] +1 == 1) {
for (int j = s[1] +1; j < e[1]; j++) {
unew[s[0]][j] = 0.25 * ( u[s[0] -1][j] + u[s[0] +1][j] + u[s[0]][j +1] + u[s[0]][j -1] - h*h*f[s[0]][j] );
}}

}


void nbxchange_and_sweep(double u[][maxn], double f[][maxn], int nx, int ny, int s, int e, double unew[][maxn], MPI_Comm comm, int nbrleft, int nbrright)
{
  MPI_Request req[4];
  MPI_Status status;
  int idx;
  double h;
  int i,j,k;

  int myid;
  MPI_Comm_rank(comm, &myid);

  h = 1.0/( (double)(nx+1) );

    /* int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, */
    /*               int source, int tag, MPI_Comm comm, MPI_Request *request); */
    /* int MPI_Isend(const void *buf, int count, MPI_Datatype datatype, int dest, */
    /* 		  int tag, MPI_Comm comm, MPI_Request *request); */

  MPI_Irecv(&u[s-1][1], ny, MPI_DOUBLE, nbrleft, 1, comm, &req[0] );
  MPI_Irecv(&u[e+1][1], ny, MPI_DOUBLE, nbrright, 2, comm, &req[1] );

  MPI_Isend(&u[e][1], ny, MPI_DOUBLE, nbrright, 1, comm, &req[2]);
  MPI_Isend(&u[s][1], ny, MPI_DOUBLE, nbrleft, 2, comm, &req[3]);
    
  /* perform purely local updates (that don't need ghosts) */
  /* 2 cols or less means all are on processor boundary */
  if( e-s+1 > 2 ){
    for(i=s+1; i<e; i++){
      for(j=1; j<ny+1; j++){
	unew[i][j] = 0.25 * ( u[i-1][j] + u[i+1][j] + u[i][j+1] + u[i][j-1]  - h*h*f[i][j] );
      }
    }
  }

  /* perform updates in j dir only for boundary cols */
  for(j=1; j<ny+1; j++){
    unew[s][j] = 0.25 * ( u[s][j+1] + u[s][j-1]  - h*h*f[s][j] );
    unew[e][j] = 0.25 * ( u[e][j+1] + u[e][j-1]  - h*h*f[e][j] );
  }

  /* int MPI_Waitany(int count, MPI_Request array_of_requests[], */
  /*      int *index, MPI_Status *status) */
  for(k=0; k < 4; k++){
    MPI_Waitany(4, req, &idx, &status);

    /* idx 0, 1 are recvs */
    switch(idx){
    case 0:
      /* printf("myid: %d case idx 0: status.MPI_TAG: %d; status.MPI_SOURCE: %d (idx: %d)\n",myid,status.MPI_TAG, status.MPI_SOURCE,idx); */
      if( nbrleft != MPI_PROC_NULL &&
	  (status.MPI_TAG != 1 || status.MPI_SOURCE != nbrleft )){
	fprintf(stderr, "Error: I don't understand the world: (tag %d; source %d)\n",
		status.MPI_TAG, status.MPI_SOURCE);
	MPI_Abort(comm, 1);
      }

      /* left ghost update completed; update local leftmost column */
      for(j=1; j<ny+1; j++){
	unew[s][j] += 0.25 * ( u[s-1][j] );
      }
      break;
    case 1:
      /* printf("myid: %d case idx 1: status.MPI_TAG: %d; status.MPI_SOURCE: %d (idx: %d)\n",myid, status.MPI_TAG, status.MPI_SOURCE,idx); */
      if(nbrright != MPI_PROC_NULL &&
	 (status.MPI_TAG != 2 || status.MPI_SOURCE != nbrright )){
	fprintf(stderr, "Error: I don't understand the world: (tag %d; source %d)\n",
		status.MPI_TAG, status.MPI_SOURCE);
	MPI_Abort(comm, 1);
      }
      /* right ghost update completed; update local rightmost
	 column */
      for(j=1; j<ny+1; j++){
	unew[e][j] += 0.25 * ( u[e+1][j] );
      }
      break;
    default:
      break;
    }
  }
  /* splitting this off to take account of case of one column assigned
     to proc -- so left and right node neighbours are ghosts so both
     the recvs must be complete*/
  for(j=1; j<ny+1; j++){
    unew[s][j] += 0.25 * ( u[s+1][j] );
    unew[e][j] += 0.25 * ( u[e-1][j] );
  }

}

double griddiff(double a[][maxn], double b[][maxn], int nx, int s, int e)
{
  double sum;
  double tmp;
  int i, j;

  sum = 0.0;

  for(i=s; i<=e; i++){
    for(j=1;j<nx+1;j++){
      tmp = (a[i][j] - b[i][j]);
      sum = sum + tmp*tmp;
    }
  }

  return sum;

}

double griddiff_2d(double a[][maxn], double b[][maxn], int nx, int s[2], int e[2]){
double sum;
double tmp;

sum = 0.0;

for(int i = s[0]; i <= e[0]; i++){
for(int j = s[1]; j <= e[1]; j++){
tmp = (a[i][j] - b[i][j]);
sum = sum + tmp*tmp;}}

return sum;

}
