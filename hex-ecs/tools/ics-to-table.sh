#!/bin/bash

# expect argument
if [ "$1" == "" ]
then
    echo "Converts line data from ics-L?-S?-Pi?.dat to an in-state/out-state table for every energy."
    echo "Usage:"
    echo "   ics-to-table <filename>"
    exit
fi

# check that the file exists
if [ ! -f "$1" ]
then
    echo "File \"$f\" does not exist."
    exit
fi

# avoid Czech decimal comma
export LC_ALL=C

# read data from the cross section file
header=$(grep "E\[Ry\]" $1)
htransitions=$(echo "$header" | cut -c16- | tr -s ' ')
vtransitions=$(echo "$htransitions" | tr -s ' ' | sed 's/ /\n/g')

# get number of transitions
Ntransitions=$(echo "$vtransitions" | wc -l)

# get number of initial states
Nistates=$(echo "$vtransitions" | cut -f1 -d- | sort | uniq | wc -l)

# get number of final states
Nfstates=$(echo "$vtransitions" | cut -f2 -d- | sort | uniq | wc -l)

# for all energies
grep -v "#" $1 | while read line
do
    # get energy and create the filename
    energy=$(echo "$line" | cut -f1 -d' ')
    filename="ics_table_E$energy.dat"
    
    # write header
    printf "#%14s" "f \\ i" > $filename
    for i in $(seq 1 $Nistates)
    do
        istate=$(echo "$htransitions" | cut -f$(( ($i - 1) * $Nfstates + 1 )) -d' ' | cut -f1 -d-)
        printf "%15s" "$istate"
    done >> $filename
    echo >> $filename
    
    # write cross sections
    for f in $(seq 1 $Nfstates)
    do
        fstate=$(echo "$htransitions" | cut -f$f -d' ' | cut -f2 -d-)
        printf "%15s" "$fstate"
        for i in $(seq 1 $Nistates)
        do
            cs=$(echo "$line" | cut -f16- | tr -s ' ' | cut -f$(( ($i - 1) * $Nfstates + $f + 1 )) -d' ')
            printf "%15s" "$cs"
        done
        echo
    done >> $filename
    
done
