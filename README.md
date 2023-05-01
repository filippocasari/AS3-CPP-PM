# ASSIGNMENT 3 particles methods (C++ version)
This assignment has been tested on MacOS.\
Running it on another OS could lead to some problems.  
## Point A
I used for plotting matplotlibcpp.
For running the program run on Unix:
```shell
bash start.sh args
```
or 
```shell
./start.sh args
```
For example:
```shell
./start.sh 300 0.1 0.1 0.01
```
For arguments help run:
```shell
./start.sh --help 
```
or
```shell
./start.sh -h
```
args are: Num particles initial temperature (double), bath temperature and verbose (-v), delta t.\
### Alternatively, you can use cmake and CMakeLists.txt. Remember to change it before running it because some libs are hardcoded. 
### Note: 
- #### For only this point initial temperature and bath temperature must be equal. 
- #### if you don't want to run it with openmp, comment pragma in the cpp file and do not include or link any libraries related to it.
- #### I changed slightly the code of matplotlibcpp to make it work. In particular I followed some suggestions given by developers in the issues section on github.



