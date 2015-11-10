#!/bin/bash

#for act in 0.0156 0.03125 0.0625 0.125 0.25 0.5 1 2 4 8 16 32 64 128
for act in 128
do
    for i in `seq 1 20`;
    do
        ../Dynamic -c 60 -u -t 0.4 -b 1.0 -a $act -i $i -s > out$act$i.out 2>&1 &
        NPROC=$(($NPROC+1))
        if [ "$NPROC" -ge 4 ]; then
            wait
            NPROC=0
        fi
    done
done

