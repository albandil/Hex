#!/bin/bash

# This test run will calculate differential cross section for the transition 1s-2p
# at impact energy Ei = 4 Ry. Altogether 60 partial waves are calculated.
# The refulting file "dwba-1s-2p.dcs" contains differential cross sections for
# angles 0 to 180 degrees, one value per one degree. This datafile has been used
# in the article
#    Benda J., Houfek K., Comput. Phys. Commun. 185 (2014) 2893-2902.
# Note that the database has been slightly changed since that time and it no
# longer implements the Born subtraction algorithm.

# remove old data
rm -f hex.db *.log *.sql *.dcs

# compute partial waves up to L = 60, use 4 processes; parallelize by 'make'
cat > Makefile <<EOF
TARGETS = 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 \
          31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60

all: \$(TARGETS)

\$(TARGETS) :: % :
	hex-dwba 1 0 2 1 4. \$@ 1000 >> L\$@.log
EOF

make -j 4

# create a new database
hex-db --new

# fill the database with the data
for SQL in *.sql; do
  hex-db --import $SQL
done

# retrieve DCS for each S and final atomic M
for Mf in -1 0 1; do
  for S in 0 1; do
    seq 0 180 | hex-db --dcs \
        --ni=1 --li=0 --mi=0 \
        --nf=2 --lf=1 --mf=$Mf \
        --S=$S --Ei=4. \
        | grep -v '#' \
        > dwba-1s-2p$Mf-S$S.dcs
  done
done

# sum all cross sections (= even columns)
paste *.dcs | awk '{
  sum=0;
  for (i=1; i<=NF; i++)
    if (i % 2 == 0)
      sum += $i;
  print $1,sum;
}' | tee dwba-1s-2p.dcs
