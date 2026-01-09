#pragma once

/* ===== Starflood Configuration ===== */

/* === Compute === */

// Uncomment to use OpenMP target constructs (if OpenMP is enabled)
#define ENABLE_OFFLOAD_SIM
#define ENABLE_OFFLOAD_VIS

// CAUTION: Experimental
// Uncomment to create a separate thread for visualization file output (requires pthreads)
// If the visualization I/O thread is not finished by the time the next frame is ready,
// the main thread will wait for it to finish. Note that this is experimental, and if, for
// some reason, the visualization I/O thread is killed or segfaults, the main thread will
// never finish waiting!
#define ENABLE_VIS_IO_THREAD

/* === General === */

// Default number of simulation timesteps to run
#define NUM_TIMESTEPS 10

// Uncomment to enable simulation timestepping (otherwise, simulation snapshots are loaded if SIM_FILENAME is defined)
// If commented, only every NUM_TIMESTEPS/OUTPUT_INTERVAL step will be processed
#define ENABLE_SIM

// Uncomment to enable visualization
#define ENABLE_VIS

#define ENABLE_TREE

#define TREE_MAX_BODIES_PER_NODE 4

//#define TREE_MAX_DEPTH
//#define TREE_MAX_NODES

// CAUTION: Experimental
// Uncomment to enable grid calculations
//#define ENABLE_GRID

// exp2(GRID_SIZE)^3
#define GRID_SIZE 8

// CAUTION: Experimental
// Uncomment to enable FFT calculations on the grid (requires fftw3)
//#define ENABLE_FFT

// CAUTION: Experimental
// Uncomment to enable FFT mulithreading (requires ENABLE_FFT and linking the appropriate fftw3 lib)
// If OpenMP is enabled, FFTW will use (at most) omp_get_max_threads() threads.
// Otherwise, FFTW will use (at most) ENABLE_FFTW_MULTITHREADING threads.
// NOTE: It appears that (with the NVIDIA HPC Compiler's libgomp, at least) fftw3_omp is
// incompatible with code that also has target (offloading) regions
#define ENABLE_FFTW_MULTITHREADING 8

// Uncomment to use an FFTW wisdomfile (FFTW will read/save OUTPUT_DIR/fftw_wisdom)
#define ENABLE_FFTW_WISDOMFILE

// CAUTION: Experimental
// Number of samples to use when distributing mass on the grid (if enabled)
#define GRID_SAMPLES 128

/* === File I/O === */

// Output directory prefix
#define OUTPUT_DIR "./out"

// Uncomment to only output files every OUTPUT_INTERVAL steps
//#define OUTPUT_INTERVAL 10

// Uncomment to load initial conditions from a snapshot
// Otherwise, initcond_generate() will be used
//#define INIT_COND_FILE "/scratch/Simulations/ICs/starflood_expanding_sphere_1.raw"

// Uncomment to write simulation snapshots to OUTPUT_DIR/SIM_FILENAME(          .raw), if ENABLE_SIM is defined
//#define SIM_FILENAME "sim/step_%04u"

// Uncomment to write visualization frames to OUTPUT_DIR/VIS_FILENAME(.pfm/.ppm/.png), if ENABLE_VIS is defined
#define VIS_FILENAME "vis/step_%04u"

// Format to use for reading/writing simulation snapshot files
// 0: RAW (no header, just a dump of sim.mem)
// 1: PLY (simple), position stored only
// 2: PLY (custom), includes custom properties
//
// CAUTION: The code for reading PLY files isn't finished yet, don't change SIM_FILE_FORMAT_READ
// from 0 (raw binary format) yet.
//
// Instead, you should convert snapshots from the raw binary format. If !enable_sim (in main.c)
// and SIM_FILE_FORMAT_READ != SIM_FILE_FORMAT_SAVE, then snapshots will be read and
// immediately saved in the format specified by SIM_FILE_FORMAT_SAVE.
#define SIM_FILE_FORMAT_READ 0
#define SIM_FILE_FORMAT_SAVE 0

// Format to use for writing visualization image files
// 0: PFM (RGB, 32-bit floating-point,     linear       )
// 1: PPM (RGB,  8-bit        integer, non-linear BT.709)
// 2: PNG (RGB,  8-bit        integer, non-linear   sRGB) (requires ENABLE_STB_IMAGE_WRITE)
#define VIS_FILE_FORMAT 2

// Uncomment to enable PNG file I/O (requires stb submodule)
// may signifigantly increase compile times
#define ENABLE_STB_IMAGE_WRITE

// Uncomment to override stb_image_write zlib compression level (higher: longer compression time, smaller files)
// 5 is the internal minimum, setting this any lower will have no additional benefit.
// Commenting this will leave the internal default (8).
#define IO_PNG_COMPRESSION_LEVEL 14

// Uncomment to log statistics to a CSV file (OUTPUT_DIR/LOG_STATISTICS)
#define LOG_STATISTICS "statistics.csv"

// Uncomment to log    sim_step() timings to a CSV file (OUTPUT_DIR/LOG_TIMINGS_SIM_STEP), if ENABLE_TIMING is defined
#define LOG_TIMINGS_SIM_STEP "timings_sim_step.csv"

// Uncomment to log visualization_draw() timings to a CSV file (OUTPUT_DIR/LOG_TIMINGS_VIS_DRAW), if ENABLE_TIMING is defined
#define LOG_TIMINGS_VIS_DRAW "timings_vis_draw.csv"

// Uncomment to log grid_fft_solve() timings to a CSV file (OUTPUT_DIR/LOG_TIMINGS_VIS_DRAW), if ENABLE_TIMING and ENABLE_FFT is defined
#define LOG_TIMINGS_FFT_EXEC "timings_fft_exec.csv"

// Uncomment to log statistics to a CSV file (OUTPUT_DIR/LOG_STATISTICS)
#define LOG_STATISTICS "statistics.csv"

/* === Logging Parameters === */

// Uncomment to enable timing
#define ENABLE_TIMING

// Uncomment to use omp_get_wtime() instead of POSIX time.h functions, requires OpenMP to be enabled
//#define TIMING_USE_OMP_GET_WTIME

// Uncomment to override the POSIX clockid (used by clock_getres(), clock_gettime(), etc.)
// When commented, it defaults to CLOCK_MONOTONIC
// Note: This has no effect when using #define TIMING_USE_OMP_GET_WTIME
//#define STARFLOOD_POSIX_CLOCKID CLOCK_REALTIME
//#define STARFLOOD_POSIX_CLOCKID CLOCK_MONOTONIC
//#define STARFLOOD_POSIX_CLOCKID CLOCK_MONOTONIC_RAW

/* === Simulation Parameters === */

// Default number of bodies in the simulation (N)
// It is recommended you set this to a power of 2
//#define NUM_BODIES 1024
//#define NUM_BODIES 2048
//#define NUM_BODIES 4096
//#define NUM_BODIES 8192
//#define NUM_BODIES 16384
//#define NUM_BODIES 32768
#define NUM_BODIES 65536
//#define NUM_BODIES 131072
//#define NUM_BODIES 262144
//#define NUM_BODIES 524288
//#define NUM_BODIES 1048576

// Uncomment to divide the solver calculations across
// multiple timesteps (faster, but less accurate)
//#define PAIRWISE_SOLVER_DECIMATION 4

// Simulation timestep size (dt)
#define TIMESTEP_SIZE 1.000e-1

// Gravitational potential softening parameter
#define EPSILON 1.000e-3

/* === Physical Constants === */

// Gravitational constant (arbitrary)
#define G 1.000

// Gravitational constant (pc * M_sol^-1 * (km/s)^2)
//#define G 4.3009172706e-3

// Gravitational constant (m^3 * kg^-1 * s^2)
// https://physics.nist.gov/cgi-bin/cuu/Value?bg
//#define G 6.6743e-11

/* === Solver === */

// Uncomment to enable smoothed-particle hydrodynamics
#define ENABLE_SPH

// Uncomment to use Kahan Summation in the solver for higher simulation accuracy
// Useful for single-precision floating-point calculations
#define ENABLE_KAHAN_SUMMATION_SOLVER

// Uncomment to use Kahan Summation in the solver for potential energy statistics
// Does not affect simulation behavior, only the accuracy of derived statistics
#define ENABLE_KAHAN_SUMMATION_ENERGY

/* === Visualization === */

// Uncomment to visualize the simulation using accumulation rasterization
// Otherwise, the gravitational potential across space is probed and plotted
#define VISUALIZATION_RASTERIZATION

// Decimate the number of bodies used for probe calculation
//#define PROBE_DECIMATION 16

// Number of samples for spatial anti-aliasing
#define SPATIAL_SAMPLES 128

// Exposure compensation, in Log2 scale (0 = the sum of all pixels per body rendered is 1.000)
#define EXPOSURE -1.000

// Orthographic projection scale, Log2 scale (more positive: more zoomed-in, more negative: more zoomed-out)
#define ORTHO_SCALE -1.000

// Uncomment to enable automatic scaling
//#define AUTOSCALING

// Uncomment to enable motion blur
//#define SHUTTER_SPEED 0.5

/* === Miscellaneous === */

// Uncomment to use floating-point double-precision
// On most platforms, single is IEEE-754 binary32 and double is IEEE-754 binary64
//#define STARFLOOD_DOUBLE_PRECISION

// Uncomment to enable memory alignment (good practice)
// 4 KiB (4*1024 = 4096 bytes) is the smallest on many architectures
// however, larger page sizes exist:
// x86 huge pages are 2 MiB (2*1024*1024 = 2097152 bytes)
// RISCV32 "megapage" is 4 MiB (4*1024*1024 = 4194304 bytes)
// There is really no point increasing this
// much before you are just wasting memory
#define STARFLOOD_ALIGNMENT 4194304

// Uncomment to override max filename length (otherwise use FILENAME_MAX)
// FILENAME_MAX is often 4096
#define STARFLOOD_FILENAME_MAX 256

// Source code information (--help/--version)
#define STARFLOOD_REPOSITORY "https://github.com/Zi7ar21/starflood"

// Version information
#define STARFLOOD_VERSION_MAJOR 0
#define STARFLOOD_VERSION_MINOR 0
#define STARFLOOD_VERSION_PATCH 1
