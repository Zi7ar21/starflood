#pragma once

/* File I/O Configuration */

// Uncomment to enable simulation dumping
//#define ENABLE_DUMPING

// Simulation dump filename pattern
#define SIMULATION_FILENAME "./out/step_%04u.data"

// Uncomment to enable visualization
#define ENABLE_VISUALIZATION

// Visualization filename pattern
#define VISUALIZATION_FILENAME "./out/step_%04u.pfm"

/* Simulation Parameters */

// Default number of bodies in the simulation (N)
//#define NUM_BODIES 1024
//#define NUM_BODIES 2048
//#define NUM_BODIES 4096
//#define NUM_BODIES 8192
//#define NUM_BODIES 16384
//#define NUM_BODIES 32768
//#define NUM_BODIES 65536
//#define NUM_BODIES 131072
#define NUM_BODIES 262144
//#define NUM_BODIES 524288
//#define NUM_BODIES 1048576

// Divide the force calculations across this many timesteps (comment to disable)
#define N_DIV 4

/*
#ifdef N_DIV
#if(0 >= N_DIV)
#error "N_DIV must be >= 1!"
#endif
#endif
*/

// Default number of timesteps to run simulation for
#define NUM_TIMESTEPS 900

// Simulation timestep size (dt)
//#define TIMESTEP_SIZE 0.03333333333333333333333333333333
//#define TIMESTEP_SIZE 0.06666666666666666666666666666667
#define TIMESTEP_SIZE 0.100

// Gravitational constant (physical units are not yet implemented)
#define G 1.0

// alternative (m^3 * kg^-1 * s^2)
// https://physics.nist.gov/cgi-bin/cuu/Value?bg
//#define G 6.6743e-11

/* Miscellaneous Configuration */

// Uncomment to enable memory alignment (good practice)
// 4 KiB (4*1024 = 4096 bytes) is the smallest on many architectures
// however, larger page sizes exist:
// x86 huge pages are 2 MiB (2*1024*1024 = 2097152 bytes)
// RISCV32 "megapage" is 4 MiB (4*1024*1024 = 4194304 bytes)
// There is really no point increasing this
// much before you are just wasting memory
#define STARFLOOD_ALIGNMENT 4194304

// Repository information (--help/--version)
#define STARFLOOD_REPOSITORY "https://github.com/Zi7ar21/starflood"

// Semantic versioning
#define STARFLOOD_VERSION_MAJOR 0
#define STARFLOOD_VERSION_MINOR 0
#define STARFLOOD_VERSION_PATCH 1
