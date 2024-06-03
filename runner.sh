#!/bin/bash

# compile the various files
gfortran -c void_parameters.f90
gfortran -c void_subroutines.f90
gfortran -c voids.f90

# link and compile 
gfortran -o voids.out void_parameters.o void_subroutines.o voids.o
./voids.out

