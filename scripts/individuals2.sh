#!/bin/bash

for i in `seq 1 100`;
do
    ./Silac -m -i $i -t 1.0 > $i.out 2>&1 &
    NPROC=$(($NPROC+1))
    if [ "$NPROC" -ge 4 ]; then
        wait
        NPROC=0
    fi
done

