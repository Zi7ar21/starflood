#pragma once

#include <stddef.h>

#include "types.h"

typedef struct {
	unsigned int N; // Number of bodies in the simulation

	unsigned int step_number; // # of steps elapsed

	void* mem; // Shared memory

	real* pot; // Potential Energy (Scalar)
	real* kin; //   Kinetic Energy (Scalar)
	real* mas; //             Mass (Scalar)
	real* pos; //       Position (3-Vector)
	real* vel; //       Velocity (3-Vector)
	real* acc; //   Acceleration (3-Vector)
} simulation_t;

// Initialize the simulation
int simulation_init(simulation_t* simulation, unsigned int N);

// Step the simulation
int simulation_step(simulation_t* simulation);

// Free the simulation
int simulation_free(simulation_t* simulation);

// Dump the simulation (write to a file)
int simulation_dump(simulation_t* simulation, const char* restrict filename);

// Load a simulation (from a file)
int simulation_load(simulation_t* simulation, const char* restrict filename);
