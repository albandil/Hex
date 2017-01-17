#!/bin/bash

# create file containing evaluation directions
rm -f tdcs.out dirs.inp
for theta2 in $(seq 1 179); do
for phi2 in $(seq 0 359); do
    echo "(45 0 1) ($theta2 $phi2 1)" >> dirs.inp
done
done

# launch two processes evaluating the TDCS
cat dirs.inp | hex-db --tdcs --ni=1 --li=0 --mi=0 --S=0 --Ei=4 > tdcs3D.singlet.out &
cat dirs.inp | hex-db --tdcs --ni=1 --li=0 --mi=0 --S=1 --Ei=4 > tdcs3D.triplet.out &
wait

# extract only interesting information from the output
grep -v "#" tdcs3D.singlet.out | awk '{ print $12, $13, $11; }' > tdcs3D.singlet.sph
grep -v "#" tdcs3D.singlet.out | awk '{ print $12, $13, $11; }' > tdcs3D.triplet.sph
paste tdcs3D.singlet.sph tdcs3D.triplet.sph | awk '{ print $1, $2, $3 + $6; }' > tdcs3D.sph

# convert into carthesian format and write to OBJ file
awk '{ print "v", $3*sin($1/57.29578)*cos($2/57.29578), 
                  $3*sin($1/57.29578)*sin($2/57.29578),
                  $3*cos($1/57.29578); }' tdcs3D.sph \
                  > tdcs3D.obj

# append quadrangle information to the OBJ file
for theta2 in $(seq 0 178); do
for phi2 in $(seq 0 359); do    
echo $theta2 $phi2
done
done | awk '{ print "f", $1*360+$2, (($1+1)%179)*360+$2, (($1+1)%179)*360+($2+1)%360, $1*360+($2+1)%360 }' >> tdcs3D.obj
