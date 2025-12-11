#pragma once

#include <stddef.h>

#include "types.h"

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
//#define NUM_BODIES 1024u
//#define NUM_BODIES 2048u
//#define NUM_BODIES 4096u
//#define NUM_BODIES 8192u
//#define NUM_BODIES 16384u
//#define NUM_BODIES 32768u
//#define NUM_BODIES 65536u
//#define NUM_BODIES 131072u
#define NUM_BODIES 262144u
//#define NUM_BODIES 524288u
//#define NUM_BODIES 1048576u

// Divide the force calculations across this many timesteps (comment to disable)
#define N_DIV 4u

/*
#ifdef N_DIV
#if(0u >= N_DIV)
#error "N_DIV must be >= 1!"
#endif
#endif
*/

// Default number of timesteps to run simulation for
#define NUM_TIMESTEPS 900u

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

// Memory alignment
// Default: 4 MiB (4 * 1024 * 1024 bytes)
#define STARFLOOD_ALIGNMENT (size_t)4194304u

// Repository information (--help/--version)
#define STARFLOOD_REPOSITORY "https://github.com/Zi7ar21/starflood"

// Semantic versioning
#define STARFLOOD_VERSION_MAJOR 0
#define STARFLOOD_VERSION_MINOR 0
#define STARFLOOD_VERSION_PATCH 1
