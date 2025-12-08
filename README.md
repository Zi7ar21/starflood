# ![Starflood](images/starflood-badge.png)

Starflood is an open-source astrophysical simulation code written in C.

## Features

- Cosmological [N-body simulation](https://en.wikipedia.org/wiki/N-body_simulation)
  - Gravity solver using pairwise summation (aka [particle-particle method](https://www.cs.cmu.edu/afs/cs/academic/class/15850c-s96/www/nbody.html#pp)), O(N²)
- Parallelization using [OpenMP](https://www.openmp.org/)
- Visualization
  - Rasterization accumulation using [atomic](https://en.wikipedia.org/wiki/Linearizability#Primitive_atomic_instructions) updates
- File I/O
  - Runs can be dumped in a binary format
  - Visualizations are saved as a frame sequence ([`.pfm` image format](https://netpbm.sourceforge.net/doc/pfm.html))
- ~~2D~~ (_New!_) **3D** _O(N*log(N))_-complexity tree gravity solver using the **[Barnes-Hut algorithm](https://en.wikipedia.org/wiki/Barnes%E2%80%93Hut_simulation)**!
- ~~2D~~ (_New!_) **3D** hydrodynamics with **[Smoothed-particle hydrodynamics (SPH)](https://en.wikipedia.org/wiki/Smoothed-particle_hydrodynamics)**!
- Simultaneous **Simulation** AND **Rendering**!
- **Profiling** using `omp_get_wtime` for `double`-precision timing (a little less granularity than a nanosecond after 6 months uptime)!

### Planned (To-do List)

- Add support for [distributed computing](https://en.wikipedia.org/wiki/Distributed_computing)
  - [Message Passing Interface (MPI)](https://en.wikipedia.org/wiki/Message_Passing_Interface)
- Add support for [Fast Multipole Method (FMM)](https://en.wikipedia.org/wiki/Fast_multipole_method) to physical solver(s)
- Add support for [GPU compute](https://en.wikipedia.org/wiki/General-purpose_computing_on_graphics_processing_units)
- Add support for [volume rendering](https://en.wikipedia.org/wiki/Volume_rendering) to visualization

## References

- _[Create Your Own Smoothed-Particle Hydrodynamics Simulation (With Python)](https://github.com/pmocz/sph-python)_ by Philip Mocz
- _[The Barnes-Hut Approximation: Efficient computation of N-body forces](https://jheer.github.io/barnes-hut/)_ by Jeffrey Heer
- _[stb: single-file public domain (or MIT licensed) libraries for C/C++](https://github.com/nothings/stb)_
  - `stb_image.h` and `stb_image_write.h` headers for image file I/O

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

[to-do]

### Building the Code

[to-do]

```sh
nvc -mp=gpu -gpu=ccnative -g -march=native -O4 -pedantic -std=c11 -Wall -Wextra -Wshadow -o starflood -lm src/main.c src/rng.c src/simulation.c src/visualization.c
```

### Mounting a `tmpfs` (Optional)

Disk I/O can be annoying during development if you aren't planning on keeping run data or large frame sqeuences around for a while.

Thankfully, Linux has an easy way to create a virtual memory filesystem (sometimes called a "[ramdisk](https://en.wikipedia.org/wiki/RAM_drive)" by Windows users) using the `mount` comand:

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

Please see [Tmpfs in The Linux Kernel documentation](https://www.kernel.org/doc/html/latest/filesystems/tmpfs.html) or `tmpfs` man page for more details.

```sh
man tmpfs
```

### Running Starflood

[to-do]

### Profiling Starflood

[to-do]
