#pragma once

#include <stddef.h>

#include "config.h"
#include "types.h"

#ifdef ENABLE_GRID
#include "grid.h"
#endif

#ifdef ENABLE_TREE
#include "tree.h"
#endif

// enumeration of simulation parameters (configuration space)
enum sim_conf {
	// Position (for reserving 3-vector components)
	_SIM_POS_0, // pos array beginning
	_SIM_POS_1, // pos array (+1*N)
	_SIM_POS_2, // pos array (+2*N)

	// Velocity (for reserving 3-vector components)
	_SIM_VEL_0, // vel array beginning
	_SIM_VEL_1, // pos array (+1*N)
	_SIM_VEL_2, // pos array (+2*N)

	// Acceleration (for reserving 3-vector components)
	_SIM_ACC_0, // acc array beginning
	_SIM_ACC_1, // pos array (+1*N)
	_SIM_ACC_2, // pos array (+2*N)

	SIM_MAS, // Mass (Scalar)

	SIM_POT, // Gravitational potential (scalar)

	SIM_KEN, // Kinetic energy (Scalar)
	SIM_PEN, // Potential energy (Scalar)

	#ifdef ENABLE_SPH
	SIM_RAD, // SPH smoothing radius (Scalar)
	SIM_RHO, // SPH density (Scalar)
	SIM_PRS, // SPH pressure (Scalar)
	#endif

	SIM_DOF, // degrees of freedom in the configuration space

	SIM_POS = _SIM_POS_0, // Position (3-Vector)
	SIM_VEL = _SIM_VEL_0, // Velocity (3-Vector)
	SIM_ACC = _SIM_ACC_0  // Acceleration (3-Vector)
};

typedef struct {
	unsigned int step_number; // Number of steps elapsed

	unsigned int N; // Number of bodies in the simulation

	void* mem; // Shared memory for the simulation

	#ifdef ENABLE_GRID
	grid_t grid;
	#endif

	#ifdef ENABLE_TREE
	tree_t tree;
	#endif
} sim_t;

// Locates a simulation parameter
real* sim_find(const sim_t* restrict sim_ptr, enum sim_conf parameter_enum);

// Initialize a simulation with N particles
int sim_init(sim_t* restrict sim_ptr, unsigned int N);

// Free a simulation
int sim_free(sim_t* restrict sim_ptr);

// Read a snapshot of the simulation state (raw binary format)
int sim_read_raw(sim_t* restrict sim_ptr, const char* restrict filename);

// Read a snapshot of the simulation state (PLY Polygon File Format)
int sim_read_ply(sim_t* restrict sim_ptr, const char* restrict filename);

// Save a snapshot of the simulation state (raw binary format)
int sim_save_raw(const sim_t* restrict sim_ptr, const char* restrict filename);

// Save a snapshot of the simulation state (PLY Polygon File Format)
int sim_save_ply(const sim_t* restrict sim_ptr, const char* restrict filename);

// Update/run the physical solvers
int sim_solv(sim_t* restrict sim_ptr);

// Update/run the simulation once (runs a timestep)
int sim_step(sim_t* restrict sim_ptr);
