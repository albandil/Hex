#!/bin/bash

rm -f log

hex-ecs --example | tee -a log

sed -i 's/b) Real knots that are exclusive to the projectile, if any. (Start from zero.)/b) Real knots that are exclusive to the projectile, if any. (Start from zero.)\n  L    0  100  101/g' example.inp

hex-ecs --input example.inp | tee -a log
