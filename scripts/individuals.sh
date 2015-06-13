#!/bin/bash

for act in 0.01 0.02 0.04 0.08 0.16 0.32 0.64 1.28 2.56 5.12 10.24;
do
    for i in `seq 1 10`;
    do
        ./me2me3 -c 60 -m -a $act -i $i > out$act_$i.txt 2>&1 &
        NPROC=$(($NPROC+1))
        if [ "$NPROC" -ge 4 ]; then
            wait
            NPROC=0
        fi
    done
done

