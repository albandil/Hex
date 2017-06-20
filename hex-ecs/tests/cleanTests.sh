#!/bin/bash

for DIR in $(ls -d */)
do

    echo "Cleaning test \"$DIR\"..."
    
    (
        cd $DIR
        rm -f *.hdf log *.log log.* *.dat *.sql *.inp *.bin *.ooc *.vtk
    )

done
