# ![Starflood](images/starflood-badge.png)

> "_Wow, it's just like the real cosmological codes!_"
> \- **Somebody** _hopefully_, (2024)

Starflood is an open-source SPH and N-body code(s) written in C/C++.

## Features

- ~~2D~~ (_New!_) **3D** _O(N*log(N))_-complexity tree gravity solver using the **[Barnes-Hut algorithm](https://en.wikipedia.org/wiki/Barnes%E2%80%93Hut_simulation)**!
- ~~2D~~ (_New!_) **3D** hydrodynamics with **[Smoothed-particle hydrodynamics (SPH)](https://en.wikipedia.org/wiki/Smoothed-particle_hydrodynamics)**!
- _Reasonably performant_ low-level C/C++ code with high-level parallelism via **[OpenMP](www.openmp.org)**!
- Simultaneous **Simulation** AND **Rendering**!
- **Profiling** using `omp_get_wtime` for `double`-precision timing (a little less granularity than a nanosecond after 6 months uptime)!

### Planned

- Combine all sources into one **clean**, **configurable**, **unified SPH/N-body** code
- **Optimize SPH neighbor search** (better than _O(N²)_)
- **Multipole Barnes-Hut** or **[Fast-Multipole Method (FMM)](https://en.wikipedia.org/wiki/Fast_multipole_method)**
- Ability to **save/load simulation data** and **render it afterwards**
- Ability to **load pre-generated initial conditions**
- Heterogenous Compute **Offloading** ~~(development is currently made difficult due to [a lack of compiler packages with build support on Arch Linux](https://bugs.archlinux.org/task/63227), and NVIDIA's HPC SDK requires AVX so that's out of the question...)~~ I now have access to an Intel Core i7-3770K, which has [AVX](https://en.wikipedia.org/wiki/Advanced_Vector_Extensions) (something which [NVIDIA's HPC SDK](https://developer.nvidia.com/hpc-sdk) unfortunately [**requires**](https://forums.developer.nvidia.com/t/nvidia-hpc-sdk-version-22-1-mpi-question/202949))
- **Volumetric**-ally **rendering SPH particles** for **eye candy visualization**
- A more **standardized build system** such as **[CMake](https://cmake.org/)**
- **[Distributed computing](https://en.wikipedia.org/wiki/Distributed_computing)** (Using **[Open MPI](https://www.open-mpi.org/)** or something)

## Credits

- Philip Mocz's Python SPH Guide
  - <https://github.com/pmocz/sph-python>
  - Used as a reference for the SPH algorithm, although the numpy syntax was a little hard to decipher...
- stb Headers
  - <https://github.com/nothings/stb>
  - `stb_image.h` and `stb_image_write.h` headers for image read/write
- Jeffrey Heer's article _The Barnes-Hut Approximation: Efficient computation of N-body forces_
  - <https://jheer.github.io/barnes-hut/>

## Running Starflood

Starflood doesn't currently have a method of installation. It's up to the user to manually choose where they want the Starflood repository and build it.

### Getting the Source Code

First, change directories to a place you like to keep your repositories. For example, I like to keep mine in `~/source`.

```sh
cd <directory you want the starflood repo in>
```

Next, clone the Starflood repository.

```sh
git clone https://github.com/Zi7ar21/starflood.git
```

Next, change directories to the newly-cloned repo.

```sh
cd starflood
```

### Configuration

Starflood is still in its infancy. As such, configuration is achieved by modifying source code. Don't be scared, I try to keep the parameters as `#define`s near the top of the source files (check `src/common.h`).

### Building from Source

Starflood uses a custom shell script for compilation. While I am flattered if you trust me, I encourage you to never run shell scripts without at least inspecting them first.

In the root of the repository there is a shell script `test.sh`, which will automatically build Starflood, clean a previous render (if one exists), and render the resulting frames with `ffmpeg`. However, you can execute each of these steps individually if you so wish. `build.sh` builds `starflood` with `g++`, and `encode.sh` renders it as a `.mp4`.

### Mounting a `tmpfs` (Optional)

Starflood renders timesteps and saves the frames to the disk, before they get rendered by `ffmpeg` as a video file. You can mount a `tmpfs` to `out` either by using the `temp.sh` script,

```sh
./temp.sh
```

or by just running `mount` directly (`2G` means **2 GiB** of storage):

```sh
sudo mount -o size=1G -t tmpfs tmpfs out
```

This will make running Starflood faster as it doesn't have to wait for images to be written to a disk. My condolences if your sysadmin does not let you mount `tmpfs`s. I plan on adding the ability to defer writing images, or writing them in a separate thread, but this hasn't been implemented yet.

### Run

The compiled binary which should be located in `<path to starflood>/build`. If you aren't using `test.sh`, you can run the binary from the root of the repo. For example, if you are in the repo directory on Linux you can just run:

```sh
./build/starflood.out
```
