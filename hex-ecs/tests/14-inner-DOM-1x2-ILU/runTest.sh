#!/bin/bash

rm -f log

hex-ecs --example | tee -a log

hex-ecs \
    --input example.inp \
    --preconditioner DOM \
    --dom-xpanels 1 \
    --dom-ypanels 2 \
    --dom-preconditioner ILU \
    --dom-sweeps 6 \
    --tolerance 1e-6 \
    | tee -a log
