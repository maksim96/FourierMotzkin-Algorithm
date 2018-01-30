#FourierMotzkin Algorithm
Tests the feasibility of a given linear program.
Finds a feasible solution, if there is any.
Tested on linux 64bit and linux bash on windows 64bit

to compile:
    g++ fourierMotzkin.cpp -o fourierMotzkin
to run:
    ./fourierMotzkin lpXY
to compile and run on all instances (assuming instances are called "lp1",...,"lp14")
    ./runAll.sh