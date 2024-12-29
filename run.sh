#!/bin/sh

for j in {1..100}
    do
    ./bubbles 7.0 0.0 10000 $j 0 &
    wait
    done
