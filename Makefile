# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (c) 2023-2026 Jacob Bingham
# You should have received a copy of the GNU General Public License along with this program.
# If not, see <https://www.gnu.org/licenses/>.

# ===== Starflood Makefile =====

# === Environment Variables ===

# AMD ROCm install directory (only if compiling for AMDGPU)
ROCM_PATH ?= /opt/rocm

# === Compiler Selection ===

# default compiler
CC ?= gcc

# override default compiler (please don't commit uncommented lines)
#CC := $(ROCM_PATH)/bin/amdclang
#CC := clang
#CC := gcc
#CC := icx
#CC := nvc

# === Optimization Flags ===

# Clang/GCC optimization flags (-O2 is the suggested default)
#CFLAGS := -O0
#CFLAGS := -Og
#CLFAGS := -O1
CFLAGS := -O2
#CFLAGS := -O3

# Clang/GCC tune for performance on compiler host machine
#  x86-64: use -march=native
# AArch64: use  -mcpu=native
#CFLAGS := $(CFLAGS) -march=native
#CFLAGS := $(CFLAGS) -mcpu=native

# Clang/GCC enable OpenMP (compiler directive-based parallelization)
CFLAGS := -fopenmp $(CFLAGS)

# GCC auto-vectorization (Clang has auto-vectorization enabled by default, GCC does not unless using -O3)
ifeq ($(CC),gcc)
#CFLAGS := -ftree-vectorize $(CFLAGS)
endif

# Clang/GCC unsafe floating-point optimizations
# Note: This makes floating-point math non-deterministic
# across different compilers/platforms/vendors
#CFLAGS := -ffast-math $(CFLAGS)
#CFLAGS := -mfpmath=sse $(CFLAGS)
CFLAGS := -fno-math-errno -fassociative-math -freciprocal-math $(CFLAGS)

# Clang/GCC link-time optimization
CFLAGS := -flto=auto $(CFLAGS)

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

# Generate debugging information (regular)
#DEBUG_CFLAGS := -g

# Clang/GCC generate debugging information with extensions for GNU Project Debugger (GDB)
DEBUG_CFLAGS := -ggdb

# NVIDIA HPC Compiler
ifeq ($(CC),nvc)
DEBUG_CFLAGS := -g

# generate debugging information without disabling optimizations
DEBUG_CFLAGS := -gopt
endif

# === Profiling-Guided Optimization ===

#CFLAGS := -fprofile-instr-generate $(CFLAGS)
#llvm-profdata merge -output=starflood.profdata default.profraw
#CFLAGS := -fprofile-instr-use=starflood.profdata

# === Regular Flags ===

# Set the C language standard
CFLAGS := $(CFLAGS) -std=c99 -Wpedantic

# Enable warnings
CFLAGS := $(CFLAGS) -Wall -Wconversion -Wextra -Wshadow

# Disable some warnings
CFLAGS := $(CFLAGS) -Wno-unused-parameter -Wno-unused-variable

# Include stb headers as if they were system headers
# (so that the compiler doesn't generate errors for stb_image_write.h)
#CFLAGS := $(CFLAGS) -isystem stb

# Link standard math library (needed for math.h functions)
LDFLAGS := -lm

# Link POSIX threads (only required for separate visualization I/O thread)
# ENABLE_VIS_IO_THREAD must also be defined in src/config.h
LDFLAGS := $(LDFLAGS) -lpthread

# Link FFTW3 (only required for grid solver)
# ENABLE_FFT must also be defined in src/config.h
#LDFLAGS := -lfftw3 $(LDFLAGS)
#LDFLAGS := -lfftw3f $(LDFLAGS)
#LDFLAGS := -lfftw3_threads -lfftw3 $(LDFLAGS)
#LDFLAGS := -lfftw3f_threads -lfftw3f $(LDFLAGS)
#LDFLAGS := -lfftw3_omp -lfftw3 $(LDFLAGS)
#LDFLAGS := -lfftw3f_omp -lfftw3f $(LDFLAGS)

# Enable debug flags
CFLAGS := $(DEBUG_CFLAGS) $(CFLAGS)

# === Make ===

BUILD_DIR := ./build

CFLAGS := $(CFLAGS) -MMD -MP -isystem stb

SOURCES := $(shell find src -name '*.c')

OBJECTS := $(SOURCES:%=$(BUILD_DIR)/%.o)

DEPENDS := $(OBJECTS:.o=.d)

.PHONY: all
all: build/starflood

.PHONY: clean
clean:
	rm -rf $(BUILD_DIR)

$(BUILD_DIR)/starflood: $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(BUILD_DIR)/%.c.o: %.c
	mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -c -o $@ $<

-include depends
