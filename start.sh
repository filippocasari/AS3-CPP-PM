#!/bin/bash -l
function cleanup {
    echo "Exiting..."
    kill "$CPP_PID"
    exit
}
# catch interrupt



# if you don't want to run it with openmp, comment pragma in the cpp file and do not include or link any libraries related to it.

# set number of openmp threads
export OMP_NUM_THREADS=2
# compile
g++ main2.cpp Particle.cpp -std=c++17 -L/Library/Frameworks/Python.framework/Versions/3.10/lib -lpython3.10 -I/usr/local/Cellar/libomp/16.0.2/include -L/usr/local/Cellar/python@3.10/3.10.10_1/Frameworks/Python.framework/Versions/3.10/include/python3.10 -I/usr/local/Cellar/python@3.10/3.10.10_1/Frameworks/Python.framework/Versions/3.10/include/python3.10 -I/usr/local/lib/python3.10/site-packages/numpy/core/include -lomp -o AS3PM
# takes args
N=$1
if [[ $1 == "--help" || $1 == "-h" ]]; then
  echo "welcome to Usage"
  echo "arguments: Number of particles, init temperature, bath temperature, verbose (-v, optional) "
  exit 0
fi
if [ "$#" -lt 3 ]; then
  echo "Insert arguments!! Info of usage: --help or -h"
  exit 1
fi
trap cleanup EXIT

INIT_TEMP=$2
THERMOSTAT_TEMP=$3
VERBOSE=$4
# print args
echo "ARG 1: $N"
echo "ARG 2: $INIT_TEMP"
echo "ARG 3: $THERMOSTAT_TEMP"
echo "ARG 4: $VERBOSE"

# check OS
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
  echo " TEST ON LINUX"

elif [[ "$OSTYPE" == "darwin"* ]]; then
  echo " TESTS ON MAC OS"
else
  echo " TESTS ON WINDOWS...check how to compile it. Linked libraries could be different"
fi
# run executable
./AS3PM "$N" "$INIT_TEMP" "$THERMOSTAT_TEMP" "$VERBOSE" &
CPP_PID=$!
# wait until program completes its job
wait $CPP_PID
echo "Particle Methods 3 Test Finished!!!"

