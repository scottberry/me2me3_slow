#!/bin/bash

for act in 1 2 4 8 16
do
    for i in `seq 1 100`;
    do
        ./Dynamic -r -u -h0.001 -t0.333 -a1.0 -b$act -s -i$i > out$act$i.out 2>&1 &
        NPROC=$(($NPROC+1))
        if [ "$NPROC" -ge 4 ]; then
            wait
            NPROC=0
        fi
    done
done


