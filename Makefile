# ===== Starflood Makefile =====

# === Environment Variables ===

# AMD ROCm install directory (only if compiling for AMDGPU)
ROCM_PATH ?= /opt/rocm

# === Compiler Selection ===

#CC := $(ROCM_PATH)/bin/amdclang
#CC := clang
#CC := gcc
#CC := icx
CC := nvc

# === Basic Optimization Flags ===

# Clang/GCC optimized for debug (default)
CFLAGS := -Og

# Clang/GCC tuned for performance on compiler host machine
#CFLAGS := -flto=auto -march=native -O2
#CFLAGS := -flto=auto -march=native -O3

# GCC flag to enable auto-vectorization (GCC doesn't enable this unless using -O3, and Clang uses its own auto-vectorization by default)
ifeq ($(CC),gcc)
#CFLAGS := -ftree-vectorize $(CFLAGS)
endif

# AMD ROCm Clang/LLVM Compiler
ifeq ($(CC),$(ROCM_PATH)/bin/amdclang)
CFLAGS := -fopenmp --offload-arch= -march=native -O2
endif

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
#CFLAGS := -gpu=ccnative,fastmath -mp=gpu -march=native -fast -O4 --diag_suppress lossy_conversion
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

# Include stb headers as if they were system headers
# (so that the compiler doesn't generate errors for stb_image_write.h)
#CFLAGS := $(CFLAGS) -isystem stb

# Link standard math library (needed for math.h functions)
LDFLAGS := -lm

# Link POSIX threads (only required for separate visualization I/O thread)
# ENABLE_VIS_IO_THREAD must also be defined in src/config.h
LDFLAGS := $(LDFLAGS) -lpthread

# Enable debug flags
CFLAGS := $(DEBUG_CFLAGS) $(CFLAGS)

# === Targets ===

OUT_DIR := ./build

SRC_DIR := ./src

SRCS := main.c initcond.c log.c simulation.c simulation_io.c solver.c visualization.c visualization_io.c

OBJS := $(addprefix $(OUT_DIR)/src/,$(addsuffix .o,$(SRCS)))

.PHONY: all
all: $(OUT_DIR)/starflood

.PHONY: clean
clean:
	rm -fdv $(OUT_DIR)/starflood $(OBJS)
	rm -fdv $(OUT_DIR)/src
#	rm -fdv $(OUT_DIR)

$(OUT_DIR)/src:
	mkdir -p $(OUT_DIR)/src

$(OUT_DIR)/starflood: $(OBJS) | $(OUT_DIR)/src
	$(CC) $(CFLAGS) -isystem stb $^ -o $@ $(LDFLAGS)

$(OUT_DIR)/src/main.c.o: $(SRC_DIR)/main.c | $(OUT_DIR)/src
	$(CC) $(CFLAGS) -c $< -o $@

#$(OUT_DIR)/src/grid.c.o: $(SRC_DIR)/grid.c | $(OUT_DIR)/src
#	$(CC) $(CFLAGS) -c $< -o $@

$(OUT_DIR)/src/initcond.c.o: $(SRC_DIR)/initcond.c | $(OUT_DIR)/src
	$(CC) $(CFLAGS) -c $< -o $@

$(OUT_DIR)/src/log.c.o: $(SRC_DIR)/log.c | $(OUT_DIR)/src
	$(CC) $(CFLAGS) -c $< -o $@

$(OUT_DIR)/src/simulation.c.o: $(SRC_DIR)/simulation.c | $(OUT_DIR)/src
	$(CC) $(CFLAGS) -c $< -o $@

$(OUT_DIR)/src/simulation_io.c.o: $(SRC_DIR)/simulation_io.c | $(OUT_DIR)/src
	$(CC) $(CFLAGS) -c $< -o $@

$(OUT_DIR)/src/solver.c.o: $(SRC_DIR)/solver.c | $(OUT_DIR)/src
	$(CC) $(CFLAGS) -c $< -o $@

$(OUT_DIR)/src/visualization.c.o: $(SRC_DIR)/visualization.c | $(OUT_DIR)/src
	$(CC) $(CFLAGS) -c $< -o $@

$(OUT_DIR)/src/visualization_io.c.o: $(SRC_DIR)/visualization_io.c | $(OUT_DIR)/src
	$(CC) $(CFLAGS) -isystem stb -c $< -o $@
