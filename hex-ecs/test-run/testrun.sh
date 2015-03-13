#!/bin/bash

# remove old data
rm -f hex.db *.dcs *.sql *.dat *.log

# create new database
hex-db --new

# compute the T-matrices using hex-ecs and insert data into the database
for L in 0 1 2 3; do
for S in 0 1; do
cat > ./ecs.inp <<EOF
# B-spline parameters 
# order        θ
  4     0.63

# real knot sequences
  L    0    0   4
  L  0.1  2.0  20
  L    3   80  78
 -1
 
# complex knot sequences
  L  80  120  41
 -1
 
# initial atomic states
  1 -1
  *
  *

# final atomic states
  1 -1
  *

# angular momenta
#  L  Pi  nL
  $L   0   4

# initial energies in Rydbergs
  L  3  3  1
 -1
  
# magnetic field
  0
EOF
hex-ecs | tee -a test-run.log
hex-db --import tmat-$L-$Pi.sql >> test-run.log
done
done

# optionally precompute some cross sections (not necessary here, DCS is being computed directly from the stored T-matrices)
hex-db --update

# extract differential cross sections using hex-db
seq 0 180 | hex-db --dcs --ni=1 --li=0 --mi=0 --nf=1 --lf=0 --mf=0 --Ei=4 --S=0 > singlet.dcs
seq 0 180 | hex-db --dcs --ni=1 --li=0 --mi=0 --nf=1 --lf=0 --mf=0 --Ei=4 --S=1 > triplet.dcs

# show graphics
gnuplot <<EOF
set title "Differential elastic ground state cross section for impact energy 4 Ry"
set xlabel "scattering angle [degrees]"
set ylabel "dcs [a0^2]"
set logscale y
plot "singlet.dcs" with lines title "singlet", \
     "triplet.dcs" with lines title "triplet", \
     "< paste singlet.dcs triplet.dcs" using 1:(\$2+\$4) with lines title "sum"
pause mouse
EOF
