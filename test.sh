#!/bin/bash

string=$1;

size=$((${#string}*2-2))

mpic++ --prefix /usr/local/share/OpenMPI -o preorder pro.cpp > /dev/null

mpirun --prefix /usr/local/share/OpenMPI -np "$size" preorder $string

rm -f preorder
