#!/usr/bin/env bash
rm a.out
gfortran gaussquad.f90
./a.out
rm -r running
mkdir running

gfortran -g Tulli13_1990.f90 progTulli12_1990.f90 -O2 ~/lapack-3.8.0/liblapack.a ~/lapack-3.8.0/librefblas.a

cp a.out running
cp fort.23 running
cp a.out running
cp parallel_script_v4.py running
cp raw_x.txt running
cp raw_w.txt running
cp job.sh running
cd running

python3 parallel_script_v4.py 



