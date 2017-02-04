#!/bin/bash

mkdir -p test-run-GPU
cd test-run-GPU

echo "Make sure that you are using a valid OpenCL device."
echo "Run \"hex-ecs --cl-list\" to get list of available platforms and devices."

/usr/bin/time --format "Walltime: %E\nCPU-seconds: %U\nMemory: %M kB" ../test-run.sh "GPU --cl-platform 0 --cl-device 0" | tee ../test-run-GPU.log

cd ..
