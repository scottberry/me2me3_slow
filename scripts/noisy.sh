#!/bin/bash

for noisy in 1 2 9 23 43 71 106 149 200 259 327 404 489 583 687 799 922 1053 1195 1346
do
    ./me2me3 -t0.333 -r -h0.001 -n $noisy -i $noisy > noisy$noisy.txt 2>&1 &
    NPROC=$(($NPROC+1))
    if [ "$NPROC" -ge 4 ]; then
        wait
        NPROC=0
    fi
done

