#!/bin/bash

rm -f log

hex-ecs --example | tee -a log

hex-ecs \
    --input example.inp \
    --preconditioner DOM \
    --dom-xpanels 2 \
    --dom-ypanels 1 \
    --dom-preconditioner ILU \
    --dom-sweeps 6 \
    --tolerance 1e-3 \
    | tee -a log
