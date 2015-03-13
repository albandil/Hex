#!/bin/bash

# remove old data
rm -f hex.db *.log *.sql *.dcs

# compute partial waves up to L = 70
for L in $(seq 0 60); do
  hex-dwba 1 0 2 1 4. $L 1000 \
    | tee -a test-run.log
done

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
}' > dwba-1s-2p.dcs
