#!/bin/bash

for G2 in 0.0 1.0 2.0 4.0 8.0
do
    ../me2me3 -r -g $G2 > out$G2.out 2>&1 &
    NPROC=$(($NPROC+1))
    if [ "$NPROC" -ge 4 ]; then
        wait
        NPROC=0
    fi
        
done

