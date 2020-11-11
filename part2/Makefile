# Makefile for Assignment 2, DTU 02635, Fall 2020
CC=gcc
CFLAGS=-Wall -std=c99 -g
LDFLAGS=-L.
LDLIBS=-lopenblas

ifeq ($(shell uname), Darwin)
	# Link against system default BLAS/LAPACK library on macOS
	LDLIBS=-llapack -lblas
endif

lssolve: matrix_io.o call_dgesv.o solve.o

.PHONY: clean
clean:
	-$(RM) *.o $(EXECUTABLE)
