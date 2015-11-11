#!/bin/bash

for i in `seq 0 80`
do
    act=$(echo "scale=6;100*e(-0.05*$i*l(10))" | bc -l)
    ./Dynamic -r -u -h0.001 -t0.333 -b1.0 -a$act -i startU$i > out$i.out 2>&1 &
    NPROC=$(($NPROC+1))
    if [ "$NPROC" -ge 4 ]; then
        wait
        NPROC=0
    fi
done

for i in `seq 0 80`
do
    act=$(echo "scale=6;100*e(-0.05*$i*l(10))" | bc -l)
    ./Dynamic -r -m -h0.001 -t0.333 -b1.0 -a$act -i startM$i > out$i.out 2>&1 &
    NPROC=$(($NPROC+1))
    if [ "$NPROC" -ge 4 ]; then
        wait
        NPROC=0
    fi
done


