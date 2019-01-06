#!/bin/bash

rm -f log*

hex-ecs --example | tee -a log

mpiexec \
    -n 2 \
    -bind-to none \
    -x OMP_NUM_THREADS=2 \
    -output-filename logs \
    hex-ecs \
        --mpi \
        --shared-scratch \
        --input example.inp \
        --preconditioner ILU \
        --lightweight-full \
    | tee -a log
