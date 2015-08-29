#!/bin/bash

for thresh in 0.2 0.4 0.6 0.8 1.0;
do
    ./Silac -m -t $thresh -r -h 0.0004 > silac$thresh.txt 2>&1 &
    NPROC=$(($NPROC+1))
    if [ "$NPROC" -ge 4 ]; then
        wait
        NPROC=0
    fi
done

for thresh in 0.2 0.4 0.6 0.8 1.0;
do
    ./Silac -t $thresh -i bal -r -h 0.0001 > bal$thresh.txt 2>&1 &
    NPROC=$(($NPROC+1))
    if [ "$NPROC" -ge 4 ]; then
        wait
        NPROC=0
    fi
done

