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

### Planned (To-do List)

- Add support for [distributed computing](https://en.wikipedia.org/wiki/Distributed_computing)
  - [Message Passing Interface (MPI)](https://en.wikipedia.org/wiki/Message_Passing_Interface)
- Add support for [GPU compute](https://en.wikipedia.org/wiki/General-purpose_computing_on_graphics_processing_units)
- Add support for [volume rendering](https://en.wikipedia.org/wiki/Volume_rendering) to visualization
- Improve physical solver(s)
  - Add support for [Barnes-Hut](https://en.wikipedia.org/wiki/Barnes–Hut_simulation) tree method
  - Add support for [Fast Multipole Method (FMM)](https://en.wikipedia.org/wiki/Fast_multipole_method)

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

#### Debug

```sh
clang -fopenmp -ggdb -march=x86-64 -Og -pedantic -std=c11 -Wall -Wconversion -Wextra -Wshadow -o starflood -lm src/*.c
```

```sh
gcc -fopenmp -ggdb -march=x86-64 -Og -pedantic -std=c11 -Wall -Wconversion -Wextra -Wshadow -o starflood -lm src/*.c
```

#### Optimized

```sh
clang -ffast-math -fopenmp -ggdb -march=native -O3 -pedantic -std=c11 -Wall -Wconversion -Wextra -Wshadow -o starflood -lm src/*.c
```

```sh
gcc -ffast-math -fopenmp -ggdb -march=native -O3 -pedantic -std=c11 -Wall -Wconversion -Wextra -Wshadow -o starflood -lm src/*.c
```

#### OpenMP Offloading

##### NVIDIA HPC SDK

```sh
nvc -gpu=ccnative -mp=gpu -fopenmp -g -march=native -O4 -pedantic -std=c11 -Wall -Wextra -Wshadow -o starflood -lm src/*.c
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

```sh
rm -f out/* && nice -9 ./starflood
```

```sh
ffmpeg -y -framerate 30 -color_primaries bt709 -colorspace bt709 -color_range full -i ./out/step_%04d.pfm -color_primaries bt709 -color_trc bt709 -colorspace bt709 -color_range limited -pix_fmt yuv420p -an -sn -c:v h264_nvenc -preset slow -tune hq -profile:v main -rc-lookahead 16383 -spatial_aq 1 -temporal_aq 1 -coder cabac -b_ref_mode middle -g 15 -bf 2 -b:v 8M -r 30 -an -sn starflood_out.mp4 && mpv --loop starflood_out.mp4
```

### Profiling Starflood

[to-do]
