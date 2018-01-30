#!/bin/bash
g++ fourierMotzkin.cpp -o fourierMotzkin
for VARIABLE in {1..14}
do
    echo "LP$VARIABLE"
	./fourierMotzkin lp$VARIABLE
    echo "===================================================="
done