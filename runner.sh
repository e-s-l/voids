#!/bin/bash

figlet -f mono9 VOIDS

# compile the various files
echo "Compiling Assest:"
gfortran -c void_parameters.f90
gfortran -c void_subroutines.f90
gfortran -c voids.f90

# link and compile
echo "Linking & Executing:"
gfortran -o voids.out void_parameters.o void_subroutines.o voids.o
./voids.out

# plot some data
echo "Plotting test data:"
gnuplot -e "file='delta.ini.dat'" data_plotter.p

echo "Process Complete."

