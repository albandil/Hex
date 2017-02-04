#!/bin/bash

mkdir -p test-run-ILU
cd test-run-ILU

/usr/bin/time --format "Walltime: %E\nCPU-seconds: %U\nMemory: %M kB" ../test-run.sh ILU | tee ../test-run-ILU.log

cd ..
