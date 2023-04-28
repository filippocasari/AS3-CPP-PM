#!/bin/bash -l
function cleanup {
    echo "Exiting..."
    kill "$PYTHON_PID"
    kill "$CPP_PID"
    exit
}
trap cleanup EXIT
export OMP_NUM_THREADS=6
export LDFLAGS="-L/opt/homebrew/opt/llvm/lib"
export CPPFLAGS="-I/opt/homebrew/opt/llvm/include"
#cmake ./CMakeLists.txt
#make

/opt/homebrew/opt/llvm/bin/clang++ $LDFLAGS $CPPFLAGS -std=c++11 -Xpreprocessor -fopenmp -lomp -lczmq -O3 -arch arm64 main.cpp Particle.cpp -o AS3PM
python3.11 plot_results.py &
PYTHON_PID=$!
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
wait $PYTHON_PID
echo "Particle Methods 3 Test Finished!!!"

