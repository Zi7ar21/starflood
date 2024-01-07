#!/bin/bash

# -Wno-missing-field-initializers -Wno-unused-parameter -Wno-unused-variable
#-mp=gpu

g++ -flto -fopenmp -march=native -g -Og -pedantic -std=c++17 -Wall -Wextra -Wshadow -Wno-missing-field-initializers -Wno-unused-parameter -Wno-unused-variable -Isrc/ -Istb/ -o build/starflood.out src/rng.cpp src/barnes-hut.cpp src/main.cpp
#nvc++ -fopenmp -march=native -g -O2 -pedantic -std=c++17 -Wall -Wextra -Wshadow  -Wno-unused-parameter -Wno-unused-variable -Isrc/ -Istb/ -o build/starflood.out src/rng.cpp -fopenmp src/barnes-hut.cpp -fopenmp src/main.cpp
#/usr/bin/time -v g++ -fopenmp -march=x86-64 -O0 -pedantic -std=c++11 -Wall -Wextra -Wshadow -Wno-missing-field-initializers -Istb/ -o build/starflood.out src/main.cpp
