#pragma once

/* ===== Starflood Configuration ===== */

/* === Compute === */

// Uncomment to utilize offloading (when using a supported toolchain with OpenMP enabled)
#define ENABLE_OFFLOADING

/* === File I/O === */

// Uncomment to only output files every N steps
#define OUTPUT_INTERVAL 5

// Comment to load SIMULATION_FILENAME instead of running simulation by default
#define ENABLE_SIMULATION

// Uncomment to enable simulation file I/O
#define SIMULATION_FILENAME "./out/sim/step_%04u.raw"

// Uncomment to enable visualization by default
#define ENABLE_VISUALIZATION

// Uncomment to enable visualization file I/O
#define VISUALIZATION_FILENAME "./out/vis/step_%04u.ppm"

/* === Logging === */

// Uncomment to enable timing
#define ENABLE_TIMING

// Uncomment to override the POSIX clockid (used by clock_getres(), clock_gettime(), etc.)
// When commented, it defaults to CLOCK_REALTIME
// Note: This has no effect when using #define TIMING_USE_OMP_GET_WTIME
#define STARFLOOD_POSIX_CLOCKID CLOCK_REALTIME
//#define STARFLOOD_POSIX_CLOCKID CLOCK_PROCESS_CPUTIME_ID
//#define STARFLOOD_POSIX_CLOCKID CLOCK_MONOTONIC_RAW

// Uncomment to use omp_get_wtime() instead of POSIX functions, requires OpenMP to be enabled
#define TIMING_USE_OMP_GET_WTIME

// Uncomment to log simulation_step() timings
#define LOG_TIMINGS_SIM_STEP

/* === Run Parameters === */

// Default number of bodies in the simulation (N)
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

// Uncomment to divide the force calculations across
// multiple timesteps (faster, but less accurate)
//#define N_DIV 4

// Default number of timesteps to run simulation for
#define NUM_TIMESTEPS 10

// Simulation timestep size (dt)
//#define TIMESTEP_SIZE 0.03333333333333333333333333333333
//#define TIMESTEP_SIZE 0.06666666666666666666666666666667
//#define TIMESTEP_SIZE 0.100
#define TIMESTEP_SIZE 1.000

// Gravitational potential softening parameter
#define EPSILON 0.001

/* === Physical Constants === */

// Gravitational constant (physical units are not yet implemented)
#define G 1.000

// alternative (pc * (km/s)^2 * M_sol^-1)
//#define G 4.3009172706e-3

// alternative (m^3 * kg^-1 * s^2)
// https://physics.nist.gov/cgi-bin/cuu/Value?bg
//#define G 6.6743e-11

/* === Visualization === */

// Number of samples for spatial anti-aliasing
#define SPATIAL_SAMPLES 512

// Exposure compensation, in Log2 scale (0 = the sum of all pixels per body rendered is 1.000)
#define EXPOSURE -3.000

// Orthographic projection scale, Log2 scale (more positive: more zoomed-in, more negative: more zoomed-out)
#define ORTHO_SCALE -1.000

// Uncomment to enable motion blur
#define SHUTTER_SPEED 0.500

/* === Miscellaneous === */

// Uncomment to enable memory alignment (good practice)
// 4 KiB (4*1024 = 4096 bytes) is the smallest on many architectures
// however, larger page sizes exist:
// x86 huge pages are 2 MiB (2*1024*1024 = 2097152 bytes)
// RISCV32 "megapage" is 4 MiB (4*1024*1024 = 4194304 bytes)
// There is really no point increasing this
// much before you are just wasting memory
#define STARFLOOD_ALIGNMENT 4194304
