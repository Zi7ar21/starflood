# ![Starflood](images/starflood-badge.png)

Starflood is an open-source astrophysical simulation code written in C.

## Features

- Cosmological [N-body simulation](https://en.wikipedia.org/wiki/N-body_simulation)
  - Gravity solver using pairwise summation (aka [particle-particle method](https://www.cs.cmu.edu/afs/cs/academic/class/15850c-s96/www/nbody.html#pp)), O(N²)
- Parallelization using [OpenMP](https://www.openmp.org/) (compiler directive-based)
  - Device/GPU compute when using a toolchain with offloading support
- Visualization
  - Accumulation rasterization using [atomic](https://en.wikipedia.org/wiki/Linearizability#Primitive_atomic_instructions) operations
  - [Spatial anti-aliasing](https://en.wikipedia.org/wiki/Spatial_anti-aliasing) ([Gaussian window function](https://en.wikipedia.org/wiki/Window_function#Gaussian_window) stochastic sampling)
- File I/O
  - Timesteps from the simulation can be saved in a raw binary format.
  - Visualization can be done during a run, or later after it has been completed.
    - The frame sequence is saved as a series of individual images, which can be encoded using tools such as [`ffmpeg`](https://ffmpeg.org/).
    - [PFM graphic image file format (`.pfm`)](https://netpbm.sourceforge.net/doc/pfm.html)

### Planned (To-do List)

- Add physically accurate initial conditions
  - [Andromeda-Milky Way collision](https://en.wikipedia.org/wiki/Andromeda–Milky_Way_collision)
  - [Molecular cloud/Stellar nursery](https://en.wikipedia.org/wiki/Molecular_cloud)
- Add support for [distributed computing](https://en.wikipedia.org/wiki/Distributed_computing)
  - [Message Passing Interface (MPI)](https://en.wikipedia.org/wiki/Message_Passing_Interface)
- Add support for [volume rendering](https://en.wikipedia.org/wiki/Volume_rendering) to visualization
- Improve physical solver(s)
  - Add support for [Barnes-Hut](https://en.wikipedia.org/wiki/Barnes–Hut_simulation) tree method
  - Add support for [Fast Multipole Method (FMM)](https://en.wikipedia.org/wiki/Fast_multipole_method)

## References

- _[Create Your Own Smoothed-Particle Hydrodynamics Simulation (With Python)](https://github.com/pmocz/sph-python)_ by Philip Mocz
  - Used as a reference implementation of [SPH](https://en.wikipedia.org/wiki/Smoothed-particle_hydrodynamics), which is nice since the Wikipedia article only contains mathematical formulae/equations!
- _[Dynamics and Astrophysics of Galaxies](https://galaxiesbook.org/)_ by Jo Bovy
  - Web-based book
- _[GADGET-4 (GAlaxies with Dark matter and Gas intEracT)](https://wwwmpa.mpa-garching.mpg.de/gadget4/)_ by Volker Springel et al.
  - Used as a reference for what a well-organized astrophysical simulation codebase should look like. Also an awesome collaboration/project that was an inspiration to me.
- _[The Barnes-Hut Approximation: Efficient computation of N-body forces](https://jheer.github.io/barnes-hut/)_ by Jeffrey Heer
  - An epic article with interactive visualizations outlining the [Barnes-Hut](https://en.wikipedia.org/wiki/Barnes–Hut_simulation) algorithm.
- _[STARFORGE: Star Formation in Gaseous Environments](https://starforge.space/)_ by Mike Grudić et al.
  - Another cool collaboration/project that inspires me.
- _[PCG, A Family of Better Random Number Generators](https://www.pcg-random.org/)_ by Melissa E. O'Neill
  - Fast/high-quality/simple hash functions used for generating the simulation intial conditions and graphical visualizations.
  - [Hash Functions for GPU Rendering](https://www.jcgt.org/published/0009/03/02/) by Mark Jarzynski and Marc Olano
    - Source of the `pcg4d()` variant of PCG.
- _[stb: single-file public domain (or MIT licensed) libraries for C/C++](https://github.com/nothings/stb)_
  - Source of the `stb_image.h` and `stb_image_write.h` C headers for image file I/O. I highly recommend anybody working on a small C/C++ project like this take a look at the available stb headers.

## Documentation

### Prerequisites

Starflood was developed on a couple of Linux machines (x86 and ARM) and is intended to be ran on UNIX-like operating systems. [Since November 2017, every system on the TOP500 list is running Linux](https://www.top500.org/statistics/details/osfam/1/). The compatibility/stability with other types of systems may vary, but feel free to open an issue or pull request if you encounter any issues.

The following sections assume you are using a \*NIX system with [Git](https://git-scm.com/), [GNU Make](https://www.gnu.org/software/make/), and a C compiler (such as [Clang](https://clang.llvm.org/) or [GCC](https://www.gnu.org/software/gcc/)) installed. In addition, [`ffmpeg`](https://ffmpeg.org/) is recommended for encoding image frame sequences as playable video files.

### Obtaining the Source Code

First, clone the Starflood repository ([`--recurse-submodules` will automatically clone and initialize any of Starflood's git submodules](https://git-scm.com/docs/gitsubmodules)):

```sh
git clone --recurse-submodules https://github.com/Zi7ar21/starflood.git
```

Then, change to the directory of the repository root:

```sh
cd starflood
```

### Configuration

Before compiling, there are a few parameters at the top of [src/config.h](src/config.h) you may want to change, such as `NUM_BODIES`:

```c
// number of bodies in the simulation (N)
#define NUM_BODIES 65536
```

### Building the Code

First, create a directory for the build. The recommended default is `build` (already in the [.gitignore](.gitignore)). You can create it now (if it doesn't already exist):

```sh
mkdir -p build
```

#### GNU Make

First, make any desired changes to the [Makefile](Makefile). Some lines are commented/uncommented near the top of the file that you might want to tweak.

For parallel compilation, use `-j` to specify the number of jobs:

```sh
make -j$(nproc) all
```

(`$(nproc)` just returns the number of available processors, it can be substituted for a number or decreased to meet memory/resource limitations).

To clean up after the build (required if any modifications have been made to the source code):

```sh
make clean
```

##### Flag Description

- `-ffast-math`: Allows replacement of standard math library functions with native instructions (i.e. `sqrt()` becomes the native x86 SSE instruction `sqrtss`). **Note**: This flag may also change associativity (order of floating-point operations), causing runs to be non-deterministic across compilers and vendors!
- `-march=native`: Tells the compiler to tune generated code for the local processor on the host machine (i.e. cache size-aware optimizations, allows the use of instruction sets such as [x86 AVX](https://en.wikipedia.org/wiki/Advanced_Vector_Extensions) or [ARM NEON](https://en.wikipedia.org/wiki/ARM_architecture_family#Advanced_SIMD_(Neon))).

#### OpenMP Offloading

Starflood supports offloading to devices (i.e. coprocessors, GPUs, etc.) using OpenMP `target` directives.

##### [Intel oneAPI DPC++/C++ Compiler](https://www.intel.com/content/www/us/en/developer/tools/oneapi/dpc-compiler-documentation.html)

```sh
source /opt/intel/oneapi/setvars.sh
```

For OpenMP offloading using the [Intel OneAPI DPC++/C++ Compiler](https://www.intel.com/content/www/us/en/developer/tools/oneapi/dpc-compiler-documentation.html):

```sh
icx -fiopenmp -fopenmp-targets=spir64 -march=native -O3 -pedantic -std=c99 -Wall -Wconversion -Wextra -Wshadow -o build/starflood src/*.c -lm
```

##### [NVIDIA HPC SDK](https://developer.nvidia.com/hpc-sdk)

> **Tip**: On Arch Linux, the NVIDIA HPC SDK has a package ([`extra/nvhpc`](https://archlinux.org/packages/extra/x86_64/nvhpc/)).

For OpenMP offloading using the [NVIDIA HPC Compilers](https://docs.nvidia.com/hpc-sdk/compilers/):

```sh
nvc -gpu=ccnative -mp=gpu -fopenmp -g -march=native -O3 -pedantic -std=c99 -Wall -Wconversion -Wextra -Wshadow -o build/starflood src/*.c -lm
```

- `-mp=gpu`: Enables compilation of OpenMP target constructs for GPU execution.
- `-gpu=ccnative`: Only generate code for visible GPUs on the compiler host machine.

### Mounting a `tmpfs` (Optional)

Disk I/O can be annoying during development if you aren't planning on keeping run data or large frame sqeuences around for a while.

Thankfully, Linux has an easy way to create a virtual memory filesystem (sometimes called a "[ramdisk](https://en.wikipedia.org/wiki/RAM_drive)" by Windows users) using the `mount` comand:

```sh
mkdir -p out
```

```sh
sudo mount -o size=4G -t tmpfs tmpfs out
```

Where `size=4G` indicates a filesystem size of 4 GiB.

Remember to unmount the tmpfs when you are finished!

```sh
sudo umount out
```

Please see [Tmpfs in The Linux Kernel documentation](https://www.kernel.org/doc/html/latest/filesystems/tmpfs.html) or `tmpfs` man page for more details.

```sh
man tmpfs
```

### Running Starflood

If file I/O is enabled, please ensure the appropriate output directories exist. The default directories are `out` for statistics/timings, `out/sim` for simulation snapshots, and `out/vis` for visualizations.

```sh
mkdir -p out/sim out/vis
```

To safely stop a run prematurely, create a file named `stop` in the output directory. If `OUTPUT_DIR` is `"./out"` (defined in [src/config.h](src/config.h)):

```sh
echo "" > out/stop
```

#### Running with Nice

[`nice`](https://en.wikipedia.org/wiki/Nice_(Unix)) is a \*NIX command that can be used to run a program with a higher/lower niceness (userspace priority).

The following command will run `./build/starflood` with a niceness value 1 higher than the shell `nice` was called from:

```sh
nice -n 1 ./build/starflood
```

A higher niceness value tells the scheduler to select a lower priority, which ensures that other processes running on your system (terminal emulators, window managers, etc.) don't get starved of resources.

There shouldn't be any problems just running `starflood` by itself, but if parallelization is enabled (using all available processors) it helps ensure you can still control the system (or do other tasks) during runs.

For the most accurate profiling however, you shouldn't use too high of a niceness value (or else you run the risk of overly frequent [context switching](https://en.wikipedia.org/wiki/Context_switch) causing high variability in execution timing).

#### Viewing the Visualization

##### Encoding an Image Frame Sequence with `ffmpeg`

Look at [`encode_ffmpeg.sh`](encode_ffmpeg.sh) for examples.

##### Playing an Image Frame Sequence with `ffplay`

Look at [`ffplay.sh`](ffplay.sh) for examples.

### Profiling Starflood

#### GNU Time

A simple way to time execution is using [GNU Time](https://www.gnu.org/software/time/):

```sh
/usr/bin/time -v ./build/starflood
```

#### Linux [`perf`](https://perfwiki.github.io/)

##### Performance Counters

You can collect and measurements from [Hardware performance counters](https://en.wikipedia.org/wiki/Hardware_performance_counter) during execution (elevated privileges required for `perf stat`):

```sh
sudo perf stat -ddd sudo -u $USER ./build/starflood
```

#### [NVIDIA Nsight Compute](https://developer.nvidia.com/nsight-compute)

[to-do]

#### [Profile-guided optimization (PGO)](https://en.wikipedia.org/wiki/Profile-guided_optimization)

Compilers such as Clang/GCC support [profile-guided optimization (PGO)](https://en.wikipedia.org/wiki/Profile-guided_optimization), a technique where profiling data collected during execution can be used to further optimize generated code for the target machine.

Useful article: <https://ddmler.github.io/compiler/2018/06/29/profile-guided-optimization.html>

##### Clang

See ["Profile Guided Optimization" in the Clang Compiler User's Manual](https://clang.llvm.org/docs/UsersManual.html#profile-guided-optimization).

##### GCC

See ["3.12 Options That Control Optimization" in the GNU Manual](https://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html).
