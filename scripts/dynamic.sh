#!/bin/bash

for act in 0.0156 0.03125 0.0625 0.125 0.25 0.5 1 2 4 8 16 32 64 128 256
do
    for i in `seq 1 10`;
    do
        ./Dynamic -r -m -h0.001 -t0.333 -b1.0 -a$act -s -i$i > out$act$i.out 2>&1 &
        NPROC=$(($NPROC+1))
        if [ "$NPROC" -ge 4 ]; then
            wait
            NPROC=0
        fi
    done
done


