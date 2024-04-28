#!/bin/bash

STB_DIR="stb" # STB directory

SOURCE_DIR="src" # Starflood source directory

BUILD_DIR="build" # Build output directory

CXX="g++" # Use GCC as the compiler
#CXX="nvc++" # Use NVIDIA's compiler

#-g -Og
CXX_FLAGS="-flto -fopenmp -march=native -Ofast -pedantic -std=c++17 -Wall -Wextra -Wshadow -Wno-missing-field-initializers -Wno-unused-parameter -Wno-unused-variable" # parameters for GCC
#CXX_FLAGS="-fopenmp -march=native -mp=gpu -g -O2 -pedantic -std=c++17 -Wall -Wextra -Wshadow  -Wno-unused-parameter -Wno-unused-variable -Isrc/ -Istb/ -o build/starflood.out src/rng.cpp -fopenmp src/barnes-hut.cpp -fopenmp src/main.cpp" # parameters for NVIDIA's HPC SDK

${CXX} ${CXX_FLAGS} -I${SOURCE_DIR}/ -I${STB_DIR}/ -o ${BUILD_DIR}/starflood.out ${SOURCE_DIR}/rng.cpp ${SOURCE_DIR}/barnes-hut.cpp ${SOURCE_DIR}/main.cpp
#/usr/bin/time -v ${CXX} ${CXX_FLAGS} -I${SOURCE_DIR}/ -I${STB_DIR}/ -o ${BUILD_DIR}/starflood.out ${SOURCE_DIR}/rng.cpp ${SOURCE_DIR}/barnes-hut.cpp ${SOURCE_DIR}/main.cpp
