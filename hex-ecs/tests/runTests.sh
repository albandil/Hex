#!/bin/bash

for DIR in $(ls -d */)
do

    printf "\033[0;32mRunning test \"$DIR\"...\033[0m\n"
    
    (
        cd $DIR
        if [ -f log ]
        then
            echo "    skipping test - log file exists"
        else
            bash runTest.sh 1>/dev/null 2>log.err
        fi
        
        if [ -f ics-L0-S0-Pi0.dat ] && [ -f ics-L0-S1-Pi0.dat ]
        then
            echo "    test succeeded"
        else
            printf "    \033[0;31mtest FAILED\033[0m\n"
        fi
    )

done
