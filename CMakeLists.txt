cmake_minimum_required(VERSION 3.10)

project(AS3PM CXX)

set(SOURCES
        main.cpp
        Particle.cpp
        )

add_executable(AS3PM ${SOURCES})

target_compile_features(AS3PM PRIVATE cxx_std_17)


target_include_directories(AS3PM PRIVATE
        /usr/local/Cellar/libomp/16.0.2/include
        /usr/local/Cellar/python@3.10/3.10.10_1/Frameworks/Python.framework/Versions/3.10/include/python3.10
        /usr/local/lib/python3.10/site-packages/numpy/core/include
        )

target_link_directories(AS3PM PRIVATE
        /usr/local/Cellar/libomp/16.0.2/include
        /Library/Frameworks/Python.framework/Versions/3.10/lib
        /usr/local/Cellar/python@3.10/3.10.10_1/Frameworks/Python.framework/Versions/3.10/include/python3.10
        )

target_link_libraries(AS3PM PRIVATE
        python3.10

        )

set_target_properties(AS3PM PROPERTIES OUTPUT_NAME "AS3PM")

target_compile_options(AS3PM PRIVATE -Wall -Wextra -pedantic)
