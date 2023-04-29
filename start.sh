#!/bin/bash -l
function cleanup {
    echo "Exiting..."
    #kill "$PYTHON_PID"
    kill "$CPP_PID"
    exit
}
trap cleanup EXIT
export OMP_NUM_THREADS=2

#cmake ./CMakeLists.txt
#make
g++ main2.cpp Particle.cpp -std=c++17 -L/Library/Frameworks/Python.framework/Versions/3.10/lib -lpython3.10  -lczmq -I/usr/local/Cellar/libomp/16.0.2/include -L/usr/local/Cellar/python@3.10/3.10.10_1/Frameworks/Python.framework/Versions/3.10/include/python3.10 -I/usr/local/Cellar/python@3.10/3.10.10_1/Frameworks/Python.framework/Versions/3.10/include/python3.10 -I/usr/local/lib/python3.10/site-packages/numpy/core/include -lomp -o AS3PM
#/opt/homebrew/opt/llvm/bin/clang++ $LDFLAGS $CPPFLAGS -std=c++11 -Xpreprocessor -fopenmp -lomp -lczmq -O3 -arch arm64 main.cpp Particle.cpp -o AS3PM
#python3.8 plot_results.py &
#PYTHON_PID=$!
N=$1
INIT_TEMP=$2
THERMOSTAT_TEMP=$3
VERBOSE=$4
echo "ARG 1: $N"
echo "ARG 2: $INIT_TEMP"
echo "ARG 3: $THERMOSTAT_TEMP"
echo "ARG 4: $VERBOSE"
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
  echo " TEST ON LINUX"

elif [[ "$OSTYPE" == "darwin"* ]]; then
  echo " TESTS ON MAC OS"
else
  echo " TESTS ON WINDOWS"
fi
./AS3PM "$N" "$INIT_TEMP" "$THERMOSTAT_TEMP" "$VERBOSE" &
CPP_PID=$!
wait $CPP_PID
echo "Particle Methods 3 Test Finished!!!"

