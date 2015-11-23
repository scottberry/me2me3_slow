#!/bin/bash

#for act in 1 2 4 8 16
#for act in 1 0.5 0.25 0.125 0.0625
for act in 8
do
    for i in `seq 1 200`;
    do
        ./Dynamic -r -u -h0.001 -t0.333 -b$act -a1.0 -i $i -s > out$act$i.out 2>&1 &
        NPROC=$(($NPROC+1))
        if [ "$NPROC" -ge 4 ]; then
            wait
            NPROC=0
        fi
    done
done

