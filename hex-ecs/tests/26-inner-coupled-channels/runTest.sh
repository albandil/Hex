#!/bin/bash

rm -f log

hex-ecs --example | tee -a log

mpiexec \
    -n 4 \
    -bind-to none \
    -x OMP_NUM_THREADS=1 \
    -output-filename log \
    hex-ecs \
        --mpi \
        --groupsize 4 \
        --shared-scratch \
        --input example.inp \
        --preconditioner coupled \
        --couple-channels \
        --lightweight-full \
        --lu mumps \
    | tee -a log
