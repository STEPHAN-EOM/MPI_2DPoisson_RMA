# @file: Makefile
# @author: Chanho Eom
# @date: 23-Apr-2023
# @brief: Makefile

CC = mpicc
CFLAGS = -g -Wall #-D_CC_OVERLAP

OBJECTS2 = main2d.o jacobi.o function.o decomp1d.o
OBJECTS3 = main2d_rmav1.o jacobi.o function.o decomp1d.o
OBJECTS4 = main2d_rmav2.o jacobi.o function.o decomp1d.o

all: main2d main2d_rmav1 main2d_rmav2

##########################################

jacobi.o: jacobi.c jacobi.h poisson1d.h
	$(CC) $(CFLAGS) -c jacobi.c -o jacobi.o

function.o: function.c function.h poisson1d.h decomp1d.h
	$(CC) $(CFLAGS) -c function.c -o function.o

decomp1d.o: decomp1d.c decomp1d.h poisson1d.h
	$(CC) $(CFLAGS) -c decomp1d.c -o decomp1d.o

main2d.o: main2d.c function.h jacobi.h poisson1d.h decomp1d.h
	$(CC) $(CFLAGS) -c main2d.c -o main2d.o

main2d_rmav1.o: main2d_rmav1.c function.h jacobi.h poisson1d.h decomp1d.h
	$(CC) $(CFLAGS) -c main2d_rmav1.c -o main2d_rmav1.o

main2d_rmav2.o: main2d_rmav2.c function.h jacobi.h poisson1d.h decomp1d.h
	$(CC) $(CFLAGS) -c main2d_rmav2.c -o main2d_rmav2.o

main2d: $(OBJECTS2)
	$(CC) $(CFLAGS) -o main2d $(OBJECTS2)

main2d_rmav1: $(OBJECTS3)
	$(CC) $(CFLAGS) -o main2d_rmav1 $(OBJECTS3)
	
main2d_rmav2: $(OBJECTS4)
	$(CC) $(CFLAGS) -o main2d_rmav2 $(OBJECTS4)


## tests

tags:
	etags *.c *.h

.PHONY: clean tags tests

clean:
	$(RM) *.o main2d main2d_rmav1 main2d_rmav2  $(TESTS) TAGS tags
