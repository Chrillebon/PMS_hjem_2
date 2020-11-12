# Makefile for Assignment 2, DTU 02635, Fall 2020
CC = gcc
CFLAGS = -Wall -std=c99 -g
LDFLAGS = -L.
LDLIBS = -lopenblas -llapack
CFILES = matrix_io.c call_dgesv.c main.c
TARGET = program

$(TARGET): $(CFILES)
	$(CC) $(CFLAGS) -o $(TARGET) $(CFILES) $(LDFLAGS) $(LDLIBS)

.PHONY: clean
clean:
	-$(RM) *.o $(EXECUTABLE)
