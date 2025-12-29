# ===== Starflood Makefile =====

CC := clang
#CC := gcc
#CC := icx
#CC := nvc

# === Basic Optimization Flags ===

# Clang/GCC optimized for debug (default)
#CFLAGS := -Og

# Clang/GCC tuned for performance on compiler host machine
CFLAGS := -flto=auto -march=native -O2
#CFLAGS := -flto=auto -march=native -O3

# Intel oneAPI DPC++/C++ Compiler tuned for performance on compiler host machine with offloading
ifeq ($(CC),icx)
CFLAGS := -fiopenmp -fopenmp-targets=spir64 -march=native -O2
endif

# NVIDIA HPC Compilers tuned for performance on compiler host machine with offloading
ifeq ($(CC),nvc)
CFLAGS := -gpu=ccnative -mp=gpu -march=native -O2 --diag_suppress lossy_conversion
#CFLAGS := -gpu=ccnative -mp=gpu -march=native -fast -O2 --diag_suppress lossy_conversion
#CFLAGS := -gpu=ccnative -mp=gpu -march=native -fast -O3 --diag_suppress lossy_conversion
#CFLAGS := -gpu=ccnative -mp=gpu -march=native -fast -O4 --diag_suppress lossy_conversion
endif

# === More Optimization Flags ===
# Uncomment any of the following special flags matching your compiler

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

# Link POSIX Threads library
# (only for visualization threaded I/O)
LDFLAGS := $(LDFLAGS) -lpthread

# === Targets ===

OUT_DIR := ./build

SRC_DIR := ./src

SRCS := main.c initcond.c log.c simulation.c solver.c visualization.c

OBJS := $(addprefix $(OUT_DIR)/,$(addsuffix .o,$(SRCS)))

.PHONY: all
all: $(OUT_DIR)/starflood

.PHONY: clean
clean:
	rm -fv $(OUT_DIR)/starflood $(OBJS)

$(OUT_DIR):
	mkdir -p $(OUT_DIR)

$(OUT_DIR)/starflood: $(OBJS) | $(OUT_DIR)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

$(OUT_DIR)/main.c.o: src/main.c
	$(CC) $(CFLAGS) -isystem stb -c $< -o $@

$(OUT_DIR)/initcond.c.o: src/initcond.c
	$(CC) $(CFLAGS) -isystem stb -c $< -o $@

$(OUT_DIR)/log.c.o: src/log.c
	$(CC) $(CFLAGS) -isystem stb -c $< -o $@

$(OUT_DIR)/simulation.c.o: src/simulation.c
	$(CC) $(CFLAGS) -isystem stb -c $< -o $@

$(OUT_DIR)/solver.c.o: src/solver.c
	$(CC) $(CFLAGS) -isystem stb -c $< -o $@

$(OUT_DIR)/visualization.c.o: src/visualization.c
	$(CC) $(CFLAGS) -isystem stb -c $< -o $@
