#pragma once

#include <stddef.h>

#include "types.h"

typedef struct {
	unsigned int step_number; // Number of steps elapsed

	unsigned int N; // Number of bodies in the simulation

	void* mem; // Shared memory for the simulation

	real* pot; // Potential Energy (Scalar)
	real* kin; //   Kinetic Energy (Scalar)
	real* mas; //             Mass (Scalar)
	real* pos; //       Position (3-Vector)
	real* vel; //       Velocity (3-Vector)
	real* acc; //   Acceleration (3-Vector)
} simulation_t;

// Initialize the simulation
int simulation_init(simulation_t* restrict simulation, unsigned int N);

// Free the simulation
int simulation_free(simulation_t* restrict simulation);

// Dump the simulation (write to a file)
int simulation_read(simulation_t* restrict simulation, const char* restrict filename);

// Read a simulation (from a file)
int simulation_save(simulation_t* restrict simulation, const char* restrict filename);

// Step the simulation
int simulation_step(simulation_t* restrict simulation);
