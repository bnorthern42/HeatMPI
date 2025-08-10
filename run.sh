#!/bin/bash
rm output.csv
make clean
make cart
mpirun -hostfile hostsfile -np 1 cartesian.out 20000 2 1