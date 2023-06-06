/**
 * @file: main2d_rmav1.c
 * @author: Chanho Eom
 * @date: 20-Apr-2023
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
int nx, ny;
int myid, nprocs;

int nbrleft, nbrright, nbrup, nbrdown;
int s[2], e[2]; 
int it;
double glob_diff;
double glob_grid_diff;
double ldiff;
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
MPI_Abort(MPI_COMM_WORLD, 1);}
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

MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &comm);
MPI_Cart_coords(comm, myid, 2, coords);
MPI_Cart_shift(comm, 0, 1, &nbrleft, &nbrright);
MPI_Cart_shift(comm, 1, 1, &nbrdown, &nbrup);

decomp_2d(nx, myid, dims, coords, s, e);

printf("(myid %d with (%d, %d)), x_axis: (%d, %d), y_axis: (%d, %d), nbrleft: %d, nbrright: %d, nbrdown: %d, nbrup: %d\n", myid, coords[0], coords[1], s[0], e[0], s[1], e[1], nbrleft, nbrright, nbrdown, nbrup);



// Question.4 => Poisson Initial Condition by C. Eom
init_basic_2d(a, b, f, nx, ny, s, e);


printf("==================== Initial Condition of grid ====================\n");
print_in_order(a, MPI_COMM_WORLD, nx);
if( nprocs == 1  ){
print_grid_to_file("grid", a,  nx, ny);
print_full_grid(a, nx);
}



t1 = MPI_Wtime();

// Question.1 => Define the windows for RMA
MPI_Win wina, winb;

// Qiestion.1 => Create the windows for RMA
MPI_Win_create(&a[s[0] -1][0], maxn * maxn * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &wina);
MPI_Win_create(&b[s[0] -1][0], maxn * maxn * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &winb);


glob_diff = 1000;
for(it=0; it<maxit; it++){

// Question1. => Implementation of the 2D ghost exchange routin: MPI_Win_fence, MPI_Put, MPU_Get
exchang_2d_RMA(a, ny, s, e, MPI_COMM_WORLD, nbrleft, nbrright, nbrdown, nbrup, wina);
sweep_2d(a, f, nx, s, e, b);
exchang_2d_RMA(b, nx, s, e, MPI_COMM_WORLD, nbrleft, nbrright, nbrdown, nbrup, winb);
sweep_2d(b, f, nx, s, e, a);


ldiff = griddiff_2d(a, b, nx, s, e);
MPI_Allreduce(&ldiff, &glob_diff, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

if(myid==0 && it%10==0){
printf("(myid %d) locdiff: %lf; glob_diff: %lf\n",myid, ldiff, glob_diff);
}

if( glob_diff < tol ){      
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


double solution[maxn][maxn];
solution_grid(solution, nx, ny);

if (myid == 0){	
printf("\n================== Show me the analytic solution ==================\n");
print_full_grid(solution, nx);}


ldiff = griddiff_2d(solution, a, nx, s, e);
MPI_Allreduce(&ldiff, &glob_grid_diff, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

if (myid == 0){
printf("\n===========================================================================\n");
printf("The global difference to the analytic solution on a grid size of %d is %f\n", nx, glob_grid_diff);
}

print_in_order(a, MPI_COMM_WORLD, nx);
if( nprocs == 1  ){
print_grid_to_file("grid", a,  nx, ny);
print_full_grid(a, nx);
}

// Question.4 => Gather the grid from the processors onto rank 0
gather_grid_2d(a, nx, ny, nprocs, myid, s, e, MPI_COMM_WORLD);
if(myid == 0){
printf("\n================== The gathered grid onto rank 0 ==================\n");
print_full_grid(a, nx);
printf("\n================== Analytic Solution ====================\n");
print_full_grid(solution, nx);
}

MPI_Win_free(&wina);
MPI_Win_free(&winb);

MPI_Finalize();
  
return 0;
}
