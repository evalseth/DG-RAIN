#!/bin/sh

gcc -std=c99 -o true_solution -I$GSL_INC true_solution_parking_lot.c -L$GSL_LIB -lm -lgsl -lgslcblas 
