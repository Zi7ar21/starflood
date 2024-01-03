# ![Starflood](images/starflood-badge.png)

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
- GPU/device offloading ~~(development is currently made difficult due to [a lack of compiler packages with build support on Arch Linux](https://bugs.archlinux.org/task/63227), and NVIDIA's HPC SDK requires AVX so that's out of the question...)~~ I now have access to an Intel Core i7-3770K, which has AVX (something which NVIDIA's HPC SDK unfortunately **requires**)
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

## Running Starflood

Starflood does not require any method of installation. It's up to the user to manually choose where they want the repository and build it.

### Getting the Source Code

First, navigate to a place you like to keep your repositories. For example, I like to keep mine in `~/source`.

```sh
cd <folder you want the starflood repo in>
```

Next, clone the Starflood repository.

```sh
git clone https://github.com/Zi7ar21/starflood.git
```

Next, enter the repository.

```sh
cd starflood
```

### Configuration

Starflood is still in its infancy, as such, all of modification happens through modifying source code. Don't fear however, as I try to keep parameters as `#define`s near the top of the source files.

### Mounting a `tmpfs` (Optional)

Starflood renders timesteps and saves the frames to the disk, before they get rendered by `ffmpeg` as a video file. Since you may be rendering a lot of frames, it is desirable to render (`2g` means a `tmpfs` with a size of _2 GiB_):

```sh
sudo mount -t tmpfs -o size=2g tmpfs out
```

This will make running Starflood faster as it doesn't wait for images to be written to the disk. My condolences if your sysadmin does not let you mount `tmpfs`s. I plan on adding the ability to defer writing images or writing them in a separate thread, but for now this is not implemented.

### Building from Source

Starflood uses a custom shell script for compilation. While I would be flattered if you trust me, I encourage you to never run shell scripts without at least inspecting them first.

In the root of the repository there is a shell script `test.sh`, which will automatically build Starflood, clean a previous render (if one exists), and render the resulting frames with `ffmpeg`. However, you can execute each of these steps individually if you so wish. `build.sh` builds `starflood` with `g++`, and `encode.sh` renders it as a `.mp4`.
