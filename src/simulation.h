#pragma once

#include <stddef.h>

#include "config.h"
#include "types.h"

// Enumeration of simulation parameters (configuration space)
enum sim_conf {
	_SIM_POS_0,
	_SIM_POS_1,
	_SIM_POS_2,

	_SIM_VEL_0,
	_SIM_VEL_1,
	_SIM_VEL_2,

	_SIM_ACC_0,
	_SIM_ACC_1,
	_SIM_ACC_2,

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
} sim_t;

// Locates a simulation parameter
real* sim_find(const sim_t* restrict sim_ptr, enum sim_conf parameter_enum);

// Initialize a simulation
int sim_init(sim_t* restrict sim_ptr, unsigned int N);

// Free a simulation
int sim_free(sim_t* restrict sim_ptr);

// Read a simulation state from a file
int sim_read(sim_t* restrict sim_ptr, const char* restrict filename);

// Save a simulation state to a file
int sim_save(const sim_t* restrict sim_ptr, const char* restrict filename);

// Updates a simulation (runs a timestep)
int sim_step(sim_t* restrict sim_ptr);
