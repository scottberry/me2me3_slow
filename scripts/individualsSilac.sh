#!/bin/bash

for i in `seq 1 1000`;
do
    ./Silac -r -m -s -i $i -t1 -h0.001 > $i.out 2>&1 &
    NPROC=$(($NPROC+1))
    if [ "$NPROC" -ge 4 ]; then
        wait
        NPROC=0
    fi
done

