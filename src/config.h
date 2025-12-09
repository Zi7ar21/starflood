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
#define NUM_BODIES 65536u

// Default number of timesteps to run simulation for
//#define NUM_TIMESTEPS 450u
#define NUM_TIMESTEPS 900u

// Simulation timestep size (dt)
#define TIMESTEP_SIZE 0.06666666666666666666666666666667

// Gravitational constant (physical units are not yet implemented)
#define G (real)1.0

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
