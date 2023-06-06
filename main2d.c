/**
 * @file: main2d.c
 * @author: Chanho Eom
 * @date: 24-Apr-2023
 * */

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

#include <mpi.h>

#include "poisson1d.h"
#include "jacobi.h"
#include "function.h"
#include "decomp1d.h"

#define maxit 10000

int main(int argc, char **argv)
{
double a[maxn][maxn], b[maxn][maxn], f[maxn][maxn];
double a1[maxn][maxn], b1[maxn][maxn], f1[maxn][maxn];
double a2[maxn][maxn], b2[maxn][maxn], f2[maxn][maxn];

int nx, ny;
int myid, nprocs;

int nbrleft, nbrright, nbrup, nbrdown;
int s[2], e[2]; 
int it;
double glob_diff, glob_diff1, glob_diff2;
double glob_grid_diff, glob_grid_diff1, glob_grid_diff2 ;
double ldiff, ldiff1, ldiff2;
double t1, t2;
double tol=1.0E-11;

MPI_Init(&argc, &argv);
MPI_Comm_rank(MPI_COMM_WORLD, &myid);
MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

MPI_Comm comm;

int dims[2] = {0, 0}; 
int periods[2] = {0, 0};
int coords[2] = {0, 0};

if( myid == 0 ){
if(argc > 2){
fprintf(stderr,"---->Usage: mpirun -np <nproc> %s <nx>\n",argv[0]);
fprintf(stderr,"---->(for this code nx=ny)\n");
MPI_Abort(MPI_COMM_WORLD, 1);}w
if( argc == 2 ){
nx = atoi(argv[1]);}
if( argc == 1 ){
nx=15;}

if( nx > maxn-2 ){
fprintf(stderr,"grid size too large\n");
exit(1);}

}


MPI_Bcast(&nx, 1, MPI_INT, 0, MPI_COMM_WORLD);
ny = nx;

printf("(myid: %d) nx = %d, ny = %d\n", myid, nx, ny);


MPI_Dims_create(nprocs, 2, dims);
if (myid == 0){
printf("The paif of processors: (x, y) = (%d, %d)\n", dims[0], dims[1]);
}

init_full_grids(a, b, f);
init_full_grids(a1, b1, f1);
init_full_grids(a2, b2, f2);

MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &comm);
MPI_Cart_coords(comm, myid, 2, coords);
MPI_Cart_shift(comm, 0, 1, &nbrleft, &nbrright);
MPI_Cart_shift(comm, 1, 1, &nbrdown, &nbrup);

decomp_2d(nx, myid, dims, coords, s, e);

printf("(myid %d with (%d, %d)), x_axis: (%d, %d), y_axis: (%d, %d), nbrleft: %d, nbrright: %d, nbrdown: %d, nbrup: %d\n", myid, coords[0], coords[1], s[0], e[0], s[1], e[1], nbrleft, nbrright, nbrdown, nbrup);



// Exercise.2 Question.2 => Poisson Initial Condition by C. Eom
init_basic_2d(a, b, f, nx, ny, s, e);
init_basic_2d(a1, b1, f1, nx, ny, s, e); 
init_basic_2d(a2, b2, f2, nx, ny, s, e);


t1 = MPI_Wtime();

MPI_Win wina, winb;
MPI_Win winc, wind;

// Exercise.3 Qiestion.1 => Create thwindows for RMA
MPI_Win_create(&a1[s[0] -1][0], maxn * maxn * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &wina);
MPI_Win_create(&b1[s[0] -1][0], maxn * maxn * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &winb);

MPI_Win_create(&a2[s[0] -1][0], maxn * maxn * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &winc);
MPI_Win_create(&b2[s[0] -1][0], maxn * maxn * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &wind);


glob_diff = 1000;
glob_diff1 = 1000;
glob_diff2 = 1000;

for(it=0; it<maxit; it++){

// Exercise.2 Question.2 => 2d ghost exchange routine with Sendrecv by C. Eom
exchang3_2d(a, ny, s, e, MPI_COMM_WORLD, nbrleft, nbrright, nbrdown, nbrup);
sweep_2d(a, f, nx, s, e, b);
exchang3_2d(b, nx, s, e, MPI_COMM_WORLD, nbrleft, nbrright, nbrdown, nbrup);
sweep_2d(b, f, nx, s, e, a);

ldiff = griddiff_2d(a, b, nx, s, e);
MPI_Allreduce(&ldiff, &glob_diff, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

if(myid==0 && it%10==0){
printf("Non_RMA] (myid %d) locdiff: %lf; glob_diff: %lf\n",myid, ldiff, glob_diff);
}

if( glob_diff < tol ){      
if(myid==0){  	
printf("iterative solve converged\n");
      }      
break;
    } 
}

MPI_Barrier(MPI_COMM_WORLD);
for(it=0; it<maxit; it++){

// Exercise.3 Question.1 => Implementation of the 2D ghost exchange routin: RMA by C. Eom
exchang_2d_RMA(a1, ny, s, e, MPI_COMM_WORLD, nbrleft, nbrright, nbrdown, nbrup, wina);
sweep_2d(a1, f1, nx, s, e, b1);
exchang_2d_RMA(b1, nx, s, e, MPI_COMM_WORLD, nbrleft, nbrright, nbrdown, nbrup, winb);
sweep_2d(b1, f1, nx, s, e, a1);

ldiff1 = griddiff_2d(a1, b1, nx, s, e);
MPI_Allreduce(&ldiff1, &glob_diff1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

if(myid==0 && it%10==0){
printf("RMA] (myid %d) locdiff: %lf; glob_diff: %lf\n",myid, ldiff1, glob_diff1);
}

if( glob_diff1 < tol){
if(myid==0){
printf("iterative solve converged\n");
       }
break;
    }           
}

MPI_Barrier(MPI_COMM_WORLD);

// Exercuse,3 Question.2 => Get all processors from MPI_COMM_WORLD
MPI_Group world_group;
MPI_Comm_group(MPI_COMM_WORLD, &world_group);

// Exercise.3 Question.2 => Create two groups for PSCW
MPI_Group source[2], dest[2];

if (nbrright != MPI_PROC_NULL) {
int ranks[1] = {nbrright};
MPI_Group_incl(world_group, 1, ranks, &source[0]);}

if (nbrleft != MPI_PROC_NULL) {
int ranks[1] = {nbrleft};
MPI_Group_incl(world_group, 1, ranks, &dest[0]);}

if (nbrup != MPI_PROC_NULL) {
int ranks[1] = {nbrup};
MPI_Group_incl(world_group, 1, ranks, &source[1]);}

if (nbrdown != MPI_PROC_NULL) {
int ranks[1] = {nbrdown};
MPI_Group_incl(world_group, 1, ranks, &dest[1]);}

for(it=0; it<maxit; it++){

// Exercise.3 Question.2 => Implementation of the 2D ghost exchange routin: RMA PSCW by C. Eom
exchang_2d_RMA_pscw(a2, ny, s, e, MPI_COMM_WORLD, nbrleft, nbrright, nbrdown, nbrup, winc, source, dest);
sweep_2d(a2, f2, nx, s, e, b2);
exchang_2d_RMA_pscw(b2, nx, s, e, MPI_COMM_WORLD, nbrleft, nbrright, nbrdown, nbrup, wind, source, dest);
sweep_2d(b2, f2, nx, s, e, a2);

ldiff2 = griddiff_2d(a2, b2, nx, s, e);
MPI_Allreduce(&ldiff2, &glob_diff2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

if(myid==0 && it%10==0){
printf("RMA PSCW] (myid %d) locdiff: %lf; glob_diff: %lf\n",myid, ldiff2, glob_diff2);
}

if( glob_diff2 < tol){
if(myid==0){
printf("iterative solve converged\n");
       }
break;
     }         
}

    
t2=MPI_Wtime();  
printf("DONE! (it: %d)\n",it);

//MPI_Barrier(MPI_COMM_WORLD);

if( myid == 0 ){    
if( it == maxit ){      
fprintf(stderr,"Failed to converge\n");    
}    
printf("Run took %lf s\n",t2-t1);  
}


// Excersie.2 Question.4 => Define Analytic solution by C.Eom
double solution[maxn][maxn];
solution_grid(solution, nx, ny);


// Exercise.2 Question.4 => Compare the results from finite difference method and analytic solution
ldiff = griddiff_2d(solution, a, nx, s, e);
MPI_Allreduce(&ldiff, &glob_grid_diff, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

ldiff1 = griddiff_2d(solution, a1, nx, s, e);
MPI_Allreduce(&ldiff1, &glob_grid_diff1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

ldiff2 = griddiff_2d(solution, a2, nx, s, e);
MPI_Allreduce(&ldiff2, &glob_grid_diff2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

if (myid == 0){
printf("\n===========================================================================\n");
printf("ver. Non_RMA] The global difference to the analytic solution on a grid size of %d is %f\n", nx, glob_grid_diff);
printf("ver. RMA] The global difference to the analytic solution on a grid size of %d is %f\n", nx, glob_grid_diff1);
printf("ver. RMA PSCW] The global difference to the analytic solution on a grid size of %d is %f\n", nx, glob_grid_diff2);
}

/*
print_in_order(a, MPI_COMM_WORLD, nx);
if( nprocs == 1  ){
print_grid_to_file("grid", a,  nx, ny);
print_full_grid(a, nx);
}
*/

// Exercise.2 Question.4 => Gather the grid from the processors onto rank 0
gather_grid_2d(a, nx, ny, nprocs, myid, s, e, MPI_COMM_WORLD);
gather_grid_2d(a1, nx, ny, nprocs, myid, s, e, MPI_COMM_WORLD);
gather_grid_2d(a2, nx, ny, nprocs, myid, s, e, MPI_COMM_WORLD);

if(myid == 0){
printf("\n================== The gathered grid onto rank 0 (ver. NON_RMA) ==================\n");
print_full_grid(a, nx);
printf("\n================== The gathered grid onto rank 0 (ver. RMA) ==================\n");
print_full_grid(a1, nx);
printf("\n================== The gathered grid onto rank 0 (ver. RMA PSCW) ==================\n");
print_full_grid(a2, nx);
}

// Exercise.3 Question.3 => Compare the results from RMA Routines and Non-RMA Routine

double local1, local2;
double global1, global2;

local1 = griddiff_2d(a, a1, nx, s, e);
MPI_Allreduce(&local1, &global1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

local2 = griddiff_2d(a, a2, nx, s, e);
MPI_Allreduce(&local2, &global2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

if (myid == 0){
printf("\n===========================================================================\n");
printf("ver. RMA] The global difference to the Non_RAM routine on a grid size of %d is %f\n", nx, global1);
printf("ver. RMA PSCW] The global difference to the Non_RMA routine on a grid size of %d is %f\n", nx, global2);
}



MPI_Win_free(&wina);
MPI_Win_free(&winb);
MPI_Win_free(&winc);
MPI_Win_free(&wind);

MPI_Finalize();
  
return 0;
}
