#!/bin/bash

for i in {1..10}  ; do
    ./tsp ulysses22.txt output_ulysses22_"$i".out
    ./tsp hk48.txt output_hk48_"$i".out
    ./tsp att48.txt output_att48_"$i".out
    ./tsp gr21.txt output_gr21_"$i".out
    ./tsp berlin52.txt output_berlin52_"$i".out
    ./tsp st70.txt output_st70_"$i".out
done
