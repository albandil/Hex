#!/bin/bash

mkdir -p test-run-KPA
cd test-run-KPA

/usr/bin/time --format "Walltime: %E\nCPU-seconds: %U\nMemory: %M kB" ../test-run.sh KPA | tee ../test-run-KPA.log

cd ..
