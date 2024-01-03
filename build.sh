#!/bin/bash

# -Wno-unused-parameter -Wno-unused-variable

g++ -fopenmp -march=native -Og -g -pedantic -std=c++17 -Wall -Wextra -Wshadow -Wno-missing-field-initializers -Wno-unused-parameter -Wno-unused-variable -Istb/ -o build/starflood.out src/main.cpp
#/usr/bin/time -v g++ -fopenmp -march=x86-64 -O0 -pedantic -std=c++11 -Wall -Wextra -Wshadow -Wno-missing-field-initializers -Istb/ -o build/starflood.out src/main.cpp
