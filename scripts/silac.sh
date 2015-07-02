#!/bin/bash

for thresh in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 1.0;
do
    ./Silac -m -t $thresh > out$thresh.txt 2>&1 &
    NPROC=$(($NPROC+1))
    if [ "$NPROC" -ge 4 ]; then
        wait
        NPROC=0
    fi
done

for thresh in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 1.0;
do
    ./Silac -t $thresh -i bal > out$thresh_bal.txt 2>&1 &
    NPROC=$(($NPROC+1))
    if [ "$NPROC" -ge 4 ]; then
        wait
        NPROC=0
    fi
done

