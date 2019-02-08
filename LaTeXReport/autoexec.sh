#!/bin/bash
# nohup sh autoexec.sh  > /dev/null 2>&1 &

gfortran random_matrix.f90 -o matrix.out
gfortran random_vector.f90 -o vector.out
./matrix.out
./vector.out
mpif90 matrix_mul_vector_parallel.f90
mpirun -np 1 a.out
mpirun -np 4 a.out
mpirun -np 8 a.out
mpirun -np 16 a.out
gfortran matrix_mul_vector_serial.f90 -o b.out
./b.out
mpif90 matrix_mul_vector_serial_with_mpi.f90 -o c.out
mpirun -np 1 c.out
gfortran walltime.for -o d.out
./d.out