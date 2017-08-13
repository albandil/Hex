#!/bin/bash

rm -f log

hex-ecs --example | tee -a log

hex-ecs \
    --input example.inp \
    --preconditioner DOM \
    --dom-xpanels 2 \
    --dom-ypanels 2 \
    --dom-preconditioner KPA \
    --dom-sweeps 6 \
    --lightweight-full \
    --tolerance 1e-3 \
    | tee -a log
