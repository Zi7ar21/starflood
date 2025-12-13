# ![Starflood](images/starflood-badge.png)

Starflood is an open-source astrophysical simulation code written in C.

## Features

- Cosmological [N-body simulation](https://en.wikipedia.org/wiki/N-body_simulation)
  - Gravity solver using pairwise summation (aka [particle-particle method](https://www.cs.cmu.edu/afs/cs/academic/class/15850c-s96/www/nbody.html#pp)), O(N²)
- Parallelization using [OpenMP](https://www.openmp.org/)
- Visualization
  - Rasterization accumulation using [atomic](https://en.wikipedia.org/wiki/Linearizability#Primitive_atomic_instructions) operations
  - [Spatial anti-aliasing](https://en.wikipedia.org/wiki/Spatial_anti-aliasing) ([Gaussian window function](https://en.wikipedia.org/wiki/Window_function#Gaussian_window) stochastic sampling)
- File I/O
  - Runs can be dumped in a binary format (machine-dependent)
  - Visualizations are saved as an RGB 32-bit floating-point (RGB32F) image frame sequence ([`.pfm` image format](https://netpbm.sourceforge.net/doc/pfm.html))

### Planned (To-do List)

- Add physically accurate initial conditions
  - [Andromeda-Milky Way collision](https://en.wikipedia.org/wiki/Andromeda–Milky_Way_collision)
  - [Molecular cloud/Stellar nursery](https://en.wikipedia.org/wiki/Molecular_cloud)
- Add support for [distributed computing](https://en.wikipedia.org/wiki/Distributed_computing)
  - [Message Passing Interface (MPI)](https://en.wikipedia.org/wiki/Message_Passing_Interface)
- Add support for [GPU compute](https://en.wikipedia.org/wiki/General-purpose_computing_on_graphics_processing_units)
- Add support for [volume rendering](https://en.wikipedia.org/wiki/Volume_rendering) to visualization
- Improve physical solver(s)
  - Add support for [Barnes-Hut](https://en.wikipedia.org/wiki/Barnes–Hut_simulation) tree method
  - Add support for [Fast Multipole Method (FMM)](https://en.wikipedia.org/wiki/Fast_multipole_method)

## References

- _[Create Your Own Smoothed-Particle Hydrodynamics Simulation (With Python)](https://github.com/pmocz/sph-python)_ by Philip Mocz
  - Used as a reference implementation of [SPH](https://en.wikipedia.org/wiki/Smoothed-particle_hydrodynamics), which is nice since the Wikipedia article only contains mathematical formulae/equations!
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

### Obtaining the Source Code

First, clone the Starflood repository:

```sh
git clone https://github.com/Zi7ar21/starflood.git
```

```sh
cd starflood
```

### Configuration

There are a few parameters at the top of [src/config.h](src/config.h) worth looking at, such as `NUM_BODIES`:

```c
// number of bodies in the simulation (N)
#define NUM_BODIES 65536
```

### Building the Code

```sh
mkdir -p build
```

```sh
make -j$(nproc) all
```

#### Debug

```sh
clang -fopenmp -ggdb -Og -pedantic -std=c99 -Wall -Wconversion -Wextra -Wshadow -o build/starflood src/*.c -lm
```

```sh
gcc -fopenmp -ggdb -Og -pedantic -std=c99 -Wall -Wconversion -Wextra -Wshadow -o build/starflood src/*.c -lm
```

- `-fopenmp`: Enables OpenMP (compiler directive-based parallelization).
- `-ggdb`: Includes debugging information/symbols ideal for [GNU Project Debugger (GDB)](https://www.gnu.org/software/gdb).
- `-pedantic -Wall -Wconversion -Wextra -Wshadow`: Enables useful warnings.
- `-lm`: Links the standard math library.

#### Optimized

```sh
clang -ffast-math -fopenmp -ggdb -march=native -O3 -pedantic -std=c99 -Wall -Wconversion -Wextra -Wshadow -o build/starflood src/*.c -lm
```

```sh
gcc -ffast-math -fopenmp -ggdb -march=native -O3 -pedantic -std=c99 -Wall -Wconversion -Wextra -Wshadow -o build/starflood src/*.c -lm
```

- `-ffast-math`: Allows replacement of standard math library functions with native instructions (i.e. `sqrt()` becomes the native x86 SSE instruction `sqrtss`). **Note**: This flag may cause runs to be non-deterministic across compilers and vendors!
- `-march=native`: Generate binaries optimized for the compiler host machine (i.e. cache optimizations and enables any supported SIMD instruction sets, such as [x86 AVX](https://en.wikipedia.org/wiki/Advanced_Vector_Extensions) or [ARM NEON](https://en.wikipedia.org/wiki/ARM_architecture_family#Advanced_SIMD_(Neon))).

#### OpenMP Offloading

##### [NVIDIA HPC SDK](https://developer.nvidia.com/hpc-sdk)

> **Tip**: On Arch Linux, the NVIDIA HPC SDK has a package ([`extra/nvhpc`](https://archlinux.org/packages/extra/x86_64/nvhpc/)).

For OpenMP offloading using the [NVIDIA HPC Compilers](https://docs.nvidia.com/hpc-sdk/compilers/):

```sh
nvc -gpu=ccnative -mp=gpu -fopenmp -g -march=native -O4 -pedantic -std=c99 -Wall -Wextra -Wshadow -o build/starflood src/*.c -lm
```

- `-mp=gpu`: Enable compilation of OpenMP target directives for the GPU.
- `-gpu=ccnative`: Generates codes only for visible GPUs on the compiler host machine.

### Mounting a `tmpfs` (Optional)

Disk I/O can be annoying during development if you aren't planning on keeping run data or large frame sqeuences around for a while.

Thankfully, Linux has an easy way to create a virtual memory filesystem (sometimes called a "[ramdisk](https://en.wikipedia.org/wiki/RAM_drive)" by Windows users) using the `mount` comand:

```sh
mkdir -p out
```

```sh
sudo mount -o size=4G -t tmpfs tmpfs out
```

Where `size=1G` indicates a filesystem size of 4 GiB. You can calculate how much storage is needed for an RGB 32-bit floating-point (RGB32F) image frame sequence using the following equation:

```math
S = 4 * 3 * W * H * N
```

where:

- `S` is the size (in bytes) of the image frame sequence
- `4` is the number of bytes/float (assuming a float is 32 bits)
- `3` is the number of channels per pixel (RGB)
- `W` is the width of an image
- `H` is the height of an image
- `N` is the number of frames

To get the number of Gibibytes (GiB) of storage needed, divide `s` by `1024 * 1024` (`1048576`).

Just remember to unmount it when you are finished!

```sh
sudo umount out
```

Please see [Tmpfs in The Linux Kernel documentation](https://www.kernel.org/doc/html/latest/filesystems/tmpfs.html) or `tmpfs` man page for more details.

```sh
man tmpfs
```

### Running Starflood

This section is still a work-in-progress!

```sh
mkdir -p out
```

#### Running with Niceness

`nice` is a tool that can be used to run a program with a higher/lower niceness value. The following command will run `./build/starflood` with a niceness value 9 higher than the shell `nice` was called from:

```sh
nice -9 ./build/starflood
```

A higher niceness means the scheduler will deprioritize the process, meaning any other programs on your system will keep running smoothly.

#### Encoding the Visualization Frame Sequence using `ffmpeg`

```sh
ffmpeg -y -r 30 -framerate 30 -i ./out/step_%04d.pfm -c:v h264_nvenc -vf "scale=in_transfer=linear:out_transfer=bt709" -b:v 8M -pix_fmt yuv420p -an -sn -dn -g 15 -bf 2 -preset slow -tune hq -profile:v main -rc-lookahead 16383 -spatial_aq 1 -temporal_aq 1 -coder cabac -b_ref_mode middle starflood_out.mp4 && mpv --loop starflood_out.mp4
```

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
