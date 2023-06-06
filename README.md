# MPI EXERCISE 3

Finite difference methods 2D Poisson Implementation via MPI with RMA routines

## How to use
### 0. Load MPI module
$module load cports openmpi

### 1. Build the project
$ make clean && make

### 2. Run Q.1 script on Chuck(** RMA ** Routine)
$ mpirun -np [Number of processor] ./main2d_rmav1 [Number of Square Grids]
$ example: mpirun -np 4 ./main2d_rmav1 13

### 3. Run Q.2 script on Chuck(** RMA PSCW ** Routine)
$ mpirun -np [Number of Processor] ./main2d_rmav2 [Number of Square Grids]
$ example: mpirun -np 4 ./main2d_rmav2 13

### 4. Change the grid size in "Poisson1d.h"

### 5. Repeat step.1 to step.3

## How to compare the results from RMA and Non-RMA
### 1. Run Q.3 script on Chuck
$ mpirun -np [Number of Processor] ./main2d [Number of Square Grids]
$ example: mpirun -np 4 ./main2d 13

### 2. Change the grid size in "Poisson1d.h"

### 3. Repeat step.1 to step.1
