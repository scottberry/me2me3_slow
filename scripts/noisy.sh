#!/bin/bash

for noisy in 1 2 4 8 16 32 64 128 256 512 1024
do
    ./me2me3 -t0.333 -r -h0.001 -n $noisy -i $noisy > noisy$noisy.txt 2>&1 &
    NPROC=$(($NPROC+1))
    if [ "$NPROC" -ge 4 ]; then
        wait
        NPROC=0
    fi
done

