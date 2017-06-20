#!/bin/bash

rm -f log

hex-ecs --example | tee -a log
hex-ecs \
    --input example.inp \
    --preconditioner ILU \
    --lightweight-full \
    | tee -a log
