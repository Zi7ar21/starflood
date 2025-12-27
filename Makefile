# ===== Starflood Makefile =====

#CC := clang # LLVM Clang
#CC := gcc # GNU Compiler Collection
#CC := icx # Intel oneAPI DPC++/C++ Compiler
CC := nvc # NVIDIA HPC Compiler

# === Basic Optimization Flags ===
# Please ensure only one of the following sets
# of flags matching your compiler are uncommented

# Clang/GCC optimized for debug (default)
#CFLAGS := -Og

# Clang/GCC tuned for performance on compiler host machine
#CFLAGS := -flto=auto -march=native -O2
#CFLAGS := -flto=auto -march=native -O3

# Intel tuned for performance on compiler host machine with offloading
#CFLAGS := -fiopenmp -fopenmp-targets=spir64 -march=native -O2

# NVIDIA HPC Compilers tuned for performance on compiler host machine with offloading
#ifeq ($(CC),nvc)
#CFLAGS := -gpu=ccnative -mp=gpu -march=native -O2 --diag_suppress lossy_conversion
#CFLAGS := -gpu=ccnative -mp=gpu -march=native -fast -O2 --diag_suppress lossy_conversion
#CFLAGS := -gpu=ccnative -mp=gpu -march=native -fast -O3 --diag_suppress lossy_conversion
CFLAGS := -gpu=ccnative -mp=gpu -march=native -fast -O4 --diag_suppress lossy_conversion
#endif

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

# Disable warnings
CFLAGS := $(CFLAGS) -Wno-unused-parameter -Wno-unused-variable

# Use debug flags
CFLAGS := $(DEBUG_CFLAGS) $(CFLAGS)

# Link the standard math library
LDFLAGS := -lm

# Link POSIX Threads library (only for visualization threaded I/O)
LDFLAGS := $(LDFLAGS) -lpthread

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
	$(CC) $(CFLAGS) -isystem stb -c -o $@ src/visualization.c
