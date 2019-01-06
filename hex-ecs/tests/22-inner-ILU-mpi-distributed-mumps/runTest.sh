#!/bin/bash

rm -f log*

hex-ecs --example | tee -a log

mpiexec \
    -n 4 \
    -bind-to none \
    -x OMP_NUM_THREADS=1 \
    -output-filename logs \
    hex-ecs \
        --mpi \
        --groupsize 2 \
        --shared-scratch \
        --input example.inp \
        --preconditioner ILU \
        --lu mumps \
    | tee -a log
