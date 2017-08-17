#!/bin/bash

rm -f log

hex-ecs --example | tee -a log

sed -i 's/b) Real knots that are exclusive to the projectile, if any. (Start from zero.)/b) Real knots that are exclusive to the projectile, if any. (Start from zero.)\n  L    0  100  101/g' example.inp

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
