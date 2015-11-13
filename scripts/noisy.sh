#!/bin/bash

for noisy in 1 9 43 106 200 327 489 687 922 1195
do
    ./me2me3 -t0.333 -r -h0.001 -n $noisy -i $noisy > noisy$noisy.txt 2>&1 &
    NPROC=$(($NPROC+1))
    if [ "$NPROC" -ge 16 ]; then
        wait
        NPROC=0
    fi
done

