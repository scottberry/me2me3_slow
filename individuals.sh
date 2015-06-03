#!/bin/bash

for i in `seq 1 100`;
do
    ./me2me3 -c 60 -a 1.0 -i $i > out$i.txt 2>&1 &
    NPROC=$(($NPROC+1))
    if [ "$NPROC" -ge 16 ]; then
        wait
        NPROC=0
    fi
done

