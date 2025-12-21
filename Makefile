# ===== Starflood Makefile =====

#CC := clang
#CC := gcc
CC := nvc

# === Basic Optimization Flags ===
# Please ensure only one of the following sets
# of flags matching your compiler are uncommented

# Clang/GCC optimized for debug (default)
#CFLAGS := -Og

# Clang/GCC tuned for performance on compiler host machine
#CFLAGS := -flto=auto -march=native -O3

# NVIDIA HPC Compilers tuned for performance on compiler host machine with device offloading
#CFLAGS := -gpu=ccnative -mp=gpu -march=native -O2 --diag_suppress lossy_conversion
#CFLAGS := -gpu=ccnative -mp=gpu -march=native -fast -O2 --diag_suppress lossy_conversion
#CFLAGS := -gpu=ccnative -mp=gpu -march=native -fast -O3 --diag_suppress lossy_conversion
CFLAGS := -gpu=ccnative -mp=gpu -march=native -fast -O4 --diag_suppress lossy_conversion

# === More Optimization Flags ===
# Uncomment any of the following special flags

# Clang/GCC enable OpenMP (compiler directive-based parallelization)
CFLAGS := -fopenmp $(CFLAGS)

# Clang/GCC unsafe floating-point optimizations
# Note: This makes floating-point math non-deterministic
# across different compilers/platforms/vendors
#CFLAGS := -ffast-math $(CFLAGS)

# Generate debugging information (regular)
DEBUG_CFLAGS := -g

# Clang/GCC generate debugging information with extensions for GNU Project Debugger (GDB)
#DEBUG_CFLAGS := -ggdb

# === Profiling Guided Optimization ===

#CFLAGS := -fprofile-instr-generate $(CFLAGS)
#llvm-profdata merge -output=starflood.profdata default.profraw
#CFLAGS := -fprofile-instr-use=starflood.profdata

# === Regular Flags ===

# Set the C language standard
CFLAGS := $(CFLAGS) -pedantic -std=c99

# Enable warnings
CFLAGS := $(CFLAGS) -Wall -Wconversion -Wextra -Wshadow

# Use debug flags
CFLAGS := $(DEBUG_CFLAGS) $(CFLAGS)

# Link the standard math library
LDFLAGS := -lm

# === Targets ===

.PHONY: all
all: build/starflood

.PHONY: clean
clean:
	rm -fv build/*.o build/starflood

build/starflood: build/initcond.o build/log.o build/main.o build/simulation.o build/solver.o build/visualization.o
	$(CC) $(CFLAGS) build/*.o -o $@ $(LDFLAGS)

build/initcond.o:
	$(CC) $(CFLAGS) -c -o $@ src/initcond.c

build/log.o:
	$(CC) $(CFLAGS) -c -o $@ src/log.c

build/main.o:
	$(CC) $(CFLAGS) -c -o $@ src/main.c

build/simulation.o:
	$(CC) $(CFLAGS) -c -o $@ src/simulation.c

build/solver.o:
	$(CC) $(CFLAGS) -c -o $@ src/solver.c

build/visualization.o:
	$(CC) $(CFLAGS) -c -o $@ src/visualization.c
