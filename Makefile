# Starflood

#CC := clang
#CC := gcc
CC := nvc

# Clang/GCC optimized for debug (default)
#CFLAGS := -Og

# Clang tuned for native hardware
#CFLAGS := -flto=auto -march=native -O2

# GCC   tuned for native hardware
#CFLAGS := -flto=auto -ftree-vectorize -march=native -O2

# NVIDIA HPC SDK tuned for native hardware with device offloading
CFLAGS := -gpu=ccnative -mp=gpu -march=native -O4 --diag_suppress lossy_conversion

# non-deterministic, but faster performance
#CFLAGS := -ffast-math $(CFLAGS)

# Enable OpenMP (compiler directive-based parallelization)
CFLAGS := -fopenmp $(CFLAGS)

# Generate debugging information (regular)
DEBUG_CFLAGS := -g 
# Generate debugging information with extensions for GNU Project Debugger (GDB)
#DEBUG_CFLAGS := -ggdb




# Set the C language standard
CFLAGS := $(CFLAGS) -pedantic -std=c99
# Enable warnings
CFLAGS := $(CFLAGS) -Wall -Wconversion -Wextra -Wshadow
# Use debug flags
CFLAGS := $(DEBUG_CFLAGS) $(CFLAGS)

# Link the standard math library
LDFLAGS := -lm

.PHONY: all
all: build/starflood

build/starflood: build/main.o build/simulation.o build/solver.o build/visualization.o
	$(CC) $(CFLAGS) build/*.o -o $@ $(LDFLAGS)

build/main.o:
	$(CC) $(CFLAGS) -c -o $@ src/main.c

build/simulation.o:
	$(CC) $(CFLAGS) -c -o $@ src/simulation.c

build/solver.o:
	$(CC) $(CFLAGS) -c -o $@ src/solver.c

build/visualization.o:
	$(CC) $(CFLAGS) -c -o $@ src/visualization.c

.PHONY: clean
clean:
	rm -fv build/*.o build/starflood
