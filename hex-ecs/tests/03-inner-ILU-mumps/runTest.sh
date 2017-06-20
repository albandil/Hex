#!/bin/bash

rm -f log

hex-ecs --example | tee -a log
mpiexec \
    -n 1 \
    -bind-to none \
    hex-ecs \
        --mpi \
        --input example.inp \
        --preconditioner ILU \
        --lu mumps \
    | tee -a log
