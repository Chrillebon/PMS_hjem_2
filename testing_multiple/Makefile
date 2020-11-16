# Makefile for Assignment 2, DTU 02635, Fall 2020
CC=gcc
CFLAGS=-Wall -std=c99 -g
LDFLAGS=-L.
LDLIBS=-lopenblas
CFILES = matrix_io.c testing_multiple.c
TARGET = program

ifeq ($(shell uname), Darwin)
	# Link against system default BLAS/LAPACK library on macOS
	LDLIBS=-llapack -lblas
endif

$(TARGET): $(CFILES)
	$(CC) $(CFLAGS) -o $(TARGET) $(CFILES) $(LDFLAGS) $(LDLIBS)

.PHONY: clean
clean:
	-$(RM) $(TARGET) *.o $(EXECUTABLE)
