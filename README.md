# Starflood

Starflood is an open-source SPH and N-body code(s) written in C/C++.

## Features

- [Smoothed-particle hydrodynamics (SPH)](https://en.wikipedia.org/wiki/Smoothed-particle_hydrodynamics)
- _O(N*log(N))_-complexity tree gravity solver using the [Barnes-Hut algorithm](https://en.wikipedia.org/wiki/Barnes%E2%80%93Hut_simulation) O(NlogN)
- Decently optimized low-level C/C++ code with high-level parallelism via [OpenMP](www.openmp.org)
- Semi-simultaneous simulation and rendering
- Performance profiling using `omp_get_wtime` for `double`-precision timing (a little less precision than a nanosecond after 6 months uptime)

### Planned

- Combine all sources into one clean, configurable, 2D/3D SPH/N-body code
- Optimize SPH neighbor search to be better than _O(N²)_
- Multipole Barnes-Hut or Fast-Multipole Method
- Ability to save/load simulation data and render it afterwards
- Ability to load pre-generated initial conditions
- GPU/device offloading (development is currently made difficult due to [a lack of compiler packages with build support on Arch Linux](https://bugs.archlinux.org/task/63227), and NVIDIA's HPC SDK requires AVX so that's out of the question...)
- Volumetrically rendering SPH particles for eye candy visualization
- A more standardized build system such as [CMake](https://cmake.org/)

## Credits

- Philip Mocz's Python SPH Guide
  -  <https://github.com/pmocz/sph-python>
  - Used as a reference for the SPH algorithm, although the numpy syntax was a little hard to decipher...
- stb https://github.com/nothings/stb
  - `stb_image.h` and `stb_image_write.h` headers for image read/write
- Jeffrey Heer's article _The Barnes-Hut Approximation: Efficient computation of N-body forces_
  - <https://jheer.github.io/barnes-hut/>

## Installation

At this point in time, Starflood doesn't have a method of installion. It's up to you to manually choose where you want the repository and build it.

### Cloning the Repository

First, navigate to the folder you like to keep your repositories in. For example, I like to keep mine in `~/source`.

```sh
cd <folder you want the starflood repo in>
```

Next, clone this repository.

```sh
git clone https://github.com/Zi7ar21/starflood.git
```

Next, enter the repository

```sh
cd starflood
```

### Building Starforge

Starforge uses CMake. Outlined below are some steps to build the project in a folder inside the repository. `<starforge repository>/build` is a part of the `.gitignore`, so no need to worry about accidentally commiting all your binaries.

First, create a folder to build Starforge in.

```sh
mkdir build
```

Then enter the newly created build directory.

```sh
cd build
```

Then, run `cmake` to configure the Starforge project.

```sh
cmake ..
```

This next part is platform-dependent, but by default on most major Linux distributions CMake will create a `Makefile` so you can build the project simply by running `make`.

```sh
make
```

Congratulations! You (should) have just successfully compiled Starflood.

### Running

Starforge can be executed simply by running the compiled executable. Again, this may vary depending on the platform but on most major Linux distributions you can simply run `Starflood`.

```sh
./Starflood
```

If this does not work, your compiler may not have set the execute bit. This can be done by running `chmod +x ./Starflood`. If your sysadmin sucks, you may need to run this as a privilaged user (`sudo chmod +x ./Starflood`, `doas chmod +x ./Starflood`, or your system's equivalent). I think I have ran into this issue once, but it really shouldn't be an issue on most major Linux distributions.

On Windows the compiled binary probably has `.exe` appended to it.

```cmd
Starflood.exe
```
